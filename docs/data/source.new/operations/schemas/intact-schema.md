---
id: schemas-intact-schema
title: "IntAct Molecular Interaction Database Schema"
category: schemas
parent: _index.md
last_updated: 2026-01-22
status: migrated
tags: [schema, protein-interaction, molecular-interaction, psi-mi, imex, embl-ebi]
---

**Parent:** [Schema Documentation](./_index.md)

# IntAct Molecular Interaction Database Schema

**Document ID:** INTACT-SCHEMA
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Source Version:** IntAct Release 2025-12

---

## TL;DR

IntAct is a freely available, curated molecular interaction database from EMBL-EBI, serving as a founding member of the IMEx (International Molecular Exchange) consortium. It contains 1.2+ million binary interactions from 122,000+ publications, with deep annotation using PSI-MI controlled vocabularies. Data is available in PSI-MI XML 2.5/3.0 and PSIMITAB formats. IntAct focuses on high-quality curation with experimental evidence codes and interaction detection methods.

---

## Database Statistics (Release 2025-12)

| Metric | Count |
|--------|-------|
| **Binary Interactions** | 1,233,546 |
| **Interactor Proteins** | 147,892 |
| **Publications** | 122,341 |
| **Organisms** | 623 |
| **Human Interactions** | ~450,000 |
| **Human-specific Proteins** | ~19,500 |

---

## License

### Data License
**License:** Creative Commons Attribution 4.0 (CC BY 4.0)
**Attribution Required:** Yes
**Commercial Use:** Allowed

### Software License
**License:** Apache License 2.0
**URL:** https://github.com/EBI-IntAct

---

## Portal Access

**Web Interface:** https://www.ebi.ac.uk/intact/
**FTP Downloads:** https://ftp.ebi.ac.uk/pub/databases/intact/current/
**GitHub:** https://github.com/EBI-IntAct

---

## Data Formats

### 1. PSI-MITAB (Primary Exchange Format)

Tab-delimited format with standardized columns based on PSI-MI specification.

#### MITAB 2.5 Columns (15 columns)

| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 1 | ID(s) interactor A | Unique identifier | uniprotkb:P04637 |
| 2 | ID(s) interactor B | Unique identifier | uniprotkb:Q00987 |
| 3 | Alt. ID(s) A | Alternative identifiers | intact:EBI-366083 |
| 4 | Alt. ID(s) B | Alternative identifiers | intact:EBI-389768 |
| 5 | Alias(es) A | Gene names, synonyms | uniprotkb:TP53(gene name) |
| 6 | Alias(es) B | Gene names, synonyms | uniprotkb:MDM2(gene name) |
| 7 | Interaction detection method | PSI-MI term | psi-mi:"MI:0018"(two hybrid) |
| 8 | Publication 1st author | Author name | Momand J (1992) |
| 9 | Publication ID | PubMed ID | pubmed:1535557 |
| 10 | Taxid A | NCBI taxonomy | taxid:9606(human) |
| 11 | Taxid B | NCBI taxonomy | taxid:9606(human) |
| 12 | Interaction type | PSI-MI term | psi-mi:"MI:0915"(physical association) |
| 13 | Source database | Database name | psi-mi:"MI:0469"(IntAct) |
| 14 | Interaction ID | Database accession | intact:EBI-77734 |
| 15 | Confidence value | Scoring | intact-miscore:0.73 |

#### MITAB 2.7 Extended Columns (42 columns)

Additional columns (16-42) include:

| Column | Name | Description |
|--------|------|-------------|
| 16 | Expansion method | Complex expansion | psi-mi:"MI:1060"(spoke) |
| 17 | Biological role A | Participant role | psi-mi:"MI:0499"(unspecified role) |
| 18 | Biological role B | Participant role | psi-mi:"MI:0499"(unspecified role) |
| 19 | Experimental role A | Experimental role | psi-mi:"MI:0496"(bait) |
| 20 | Experimental role B | Experimental role | psi-mi:"MI:0498"(prey) |
| 21 | Type interactor A | Molecule type | psi-mi:"MI:0326"(protein) |
| 22 | Type interactor B | Molecule type | psi-mi:"MI:0326"(protein) |
| 23 | Xref A | Cross-references | refseq:NP_000537 |
| 24 | Xref B | Cross-references | refseq:NP_002383 |
| 25 | Xref interaction | Interaction xrefs | go:"GO:0005515"(protein binding) |
| 26 | Annotation A | Participant annotations | comment:"phosphorylated" |
| 27 | Annotation B | Participant annotations | - |
| 28 | Annotation interaction | Interaction annotations | figure legend:"Fig. 1" |
| 29 | Host organism | Experimental host | taxid:9606(human) |
| 30 | Parameters interaction | Kinetic parameters | kd:2.0e-9 |
| 31 | Creation date | Record creation | 2005/03/17 |
| 32 | Update date | Last modification | 2024/11/20 |
| 33 | Checksum A | Interactor checksum | rogid:UcdngwpTSS6hG... |
| 34 | Checksum B | Interactor checksum | rogid:K+H1TqN7A... |
| 35 | Checksum interaction | Interaction checksum | rigid:4hL... |
| 36 | Negative | Negative interaction | false |
| 37 | Feature(s) A | Binding regions | binding region:23-56 |
| 38 | Feature(s) B | Binding regions | sufficient binding region:1-109 |
| 39 | Stoichiometry A | Stoichiometry | stoichiometry:2 |
| 40 | Stoichiometry B | Stoichiometry | stoichiometry:1 |
| 41 | Identification method A | Identification | psi-mi:"MI:0102"(mass spectrometry) |
| 42 | Identification method B | Identification | psi-mi:"MI:0102"(mass spectrometry) |

---

### 2. PSI-MI XML Format

Structured XML format for complete interaction records.

#### Sample PSI-MI XML 2.5

```xml
<?xml version="1.0" encoding="UTF-8"?>
<entrySet xmlns="http://psi.hupo.org/mi/mif"
          xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
          xsi:schemaLocation="http://psi.hupo.org/mi/mif http://psidev.info/mif254"
          level="2" version="5">
  <entry>
    <source>
      <names>
        <shortLabel>IntAct</shortLabel>
        <fullName>European Bioinformatics Institute</fullName>
      </names>
      <bibref>
        <xref>
          <primaryRef db="pubmed" dbAc="MI:0446" id="14681455"/>
        </xref>
      </bibref>
      <xref>
        <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0469"/>
      </xref>
    </source>

    <experimentList>
      <experimentDescription id="EBI-77734-1">
        <names>
          <shortLabel>momand-1992-1</shortLabel>
          <fullName>MDM2 binds p53</fullName>
        </names>
        <bibref>
          <xref>
            <primaryRef db="pubmed" dbAc="MI:0446" id="1535557"/>
          </xref>
        </bibref>
        <hostOrganismList>
          <hostOrganism ncbiTaxId="9606">
            <names>
              <shortLabel>human</shortLabel>
              <fullName>Homo sapiens</fullName>
            </names>
          </hostOrganism>
        </hostOrganismList>
        <interactionDetectionMethod>
          <names>
            <shortLabel>two hybrid</shortLabel>
            <fullName>two hybrid</fullName>
          </names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0018"/>
          </xref>
        </interactionDetectionMethod>
      </experimentDescription>
    </experimentList>

    <interactorList>
      <interactor id="EBI-366083">
        <names>
          <shortLabel>p53_human</shortLabel>
          <fullName>Cellular tumor antigen p53</fullName>
          <alias type="gene name">TP53</alias>
        </names>
        <xref>
          <primaryRef db="uniprotkb" dbAc="MI:0486" id="P04637"/>
          <secondaryRef db="refseq" dbAc="MI:0481" id="NP_000537"/>
          <secondaryRef db="ensembl" dbAc="MI:0476" id="ENSG00000141510"/>
        </xref>
        <interactorType>
          <names>
            <shortLabel>protein</shortLabel>
          </names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0326"/>
          </xref>
        </interactorType>
        <organism ncbiTaxId="9606">
          <names>
            <shortLabel>human</shortLabel>
            <fullName>Homo sapiens</fullName>
          </names>
        </organism>
        <sequence>MEEPQSDPSV...</sequence>
      </interactor>

      <interactor id="EBI-389768">
        <names>
          <shortLabel>mdm2_human</shortLabel>
          <fullName>E3 ubiquitin-protein ligase Mdm2</fullName>
          <alias type="gene name">MDM2</alias>
        </names>
        <xref>
          <primaryRef db="uniprotkb" dbAc="MI:0486" id="Q00987"/>
        </xref>
        <interactorType>
          <names>
            <shortLabel>protein</shortLabel>
          </names>
          <xref>
            <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0326"/>
          </xref>
        </interactorType>
        <organism ncbiTaxId="9606">
          <names>
            <shortLabel>human</shortLabel>
          </names>
        </organism>
      </interactor>
    </interactorList>

    <interactionList>
      <interaction id="EBI-77734">
        <names>
          <shortLabel>p53-mdm2-1</shortLabel>
        </names>
        <xref>
          <primaryRef db="intact" dbAc="MI:0469" id="EBI-77734"/>
          <secondaryRef db="imex" dbAc="MI:0670" id="IM-12345"/>
        </xref>
        <experimentList>
          <experimentRef>EBI-77734-1</experimentRef>
        </experimentList>
        <participantList>
          <participant id="EBI-77734-p1">
            <interactorRef>EBI-366083</interactorRef>
            <biologicalRole>
              <names><shortLabel>unspecified role</shortLabel></names>
              <xref><primaryRef db="psi-mi" id="MI:0499"/></xref>
            </biologicalRole>
            <experimentalRoleList>
              <experimentalRole>
                <names><shortLabel>bait</shortLabel></names>
                <xref><primaryRef db="psi-mi" id="MI:0496"/></xref>
              </experimentalRole>
            </experimentalRoleList>
          </participant>
          <participant id="EBI-77734-p2">
            <interactorRef>EBI-389768</interactorRef>
            <biologicalRole>
              <names><shortLabel>unspecified role</shortLabel></names>
              <xref><primaryRef db="psi-mi" id="MI:0499"/></xref>
            </biologicalRole>
            <experimentalRoleList>
              <experimentalRole>
                <names><shortLabel>prey</shortLabel></names>
                <xref><primaryRef db="psi-mi" id="MI:0498"/></xref>
              </experimentalRole>
            </experimentalRoleList>
          </participant>
        </participantList>
        <interactionType>
          <names><shortLabel>physical association</shortLabel></names>
          <xref><primaryRef db="psi-mi" id="MI:0915"/></xref>
        </interactionType>
      </interaction>
    </interactionList>
  </entry>
</entrySet>
```

---

## REST API (PSICQUIC)

IntAct data is accessible via the PSICQUIC (Proteomics Standard Initiative Common Query Interface) web service.

### Base URL
```
https://www.ebi.ac.uk/intact/ws/interactor/
```

### PSICQUIC REST Endpoints

| Endpoint | Description |
|----------|-------------|
| `/query/{query}` | Search interactions by MIQL query |
| `/interactor/{id}` | Get interactions for specific interactor |
| `/interaction/{id}` | Get specific interaction details |

### MIQL (Molecular Interaction Query Language)

#### Query Fields

| Field | Description | Example |
|-------|-------------|---------|
| `id` | Any identifier | id:P04637 |
| `identifier` | Interactor ID | identifier:uniprotkb:P04637 |
| `alias` | Gene/protein alias | alias:TP53 |
| `pubid` | Publication ID | pubid:1535557 |
| `pubauth` | Publication author | pubauth:Momand |
| `taxid` | NCBI taxonomy | taxidA:9606 AND taxidB:9606 |
| `species` | Species name | species:human |
| `type` | Interaction type | type:"MI:0915" |
| `detmethod` | Detection method | detmethod:"MI:0018" |
| `source` | Source database | source:intact |
| `interaction_id` | Interaction accession | interaction_id:EBI-77734 |

#### Sample MIQL Queries

```
# Find all TP53 interactions
identifier:P04637

# Human-human interactions only
taxidA:9606 AND taxidB:9606

# Two-hybrid detected interactions
detmethod:"MI:0018"

# Specific publication
pubid:1535557

# Combined query
identifier:P04637 AND taxidA:9606 AND detmethod:"MI:0018"
```

### PSICQUIC Response Formats

| Format | Content-Type | Description |
|--------|--------------|-------------|
| `tab25` | text/plain | MITAB 2.5 (15 columns) |
| `tab27` | text/plain | MITAB 2.7 (42 columns) |
| `xml25` | application/xml | PSI-MI XML 2.5 |
| `count` | text/plain | Result count only |
| `json` | application/json | JSON format |

---

## Field Dictionary

### PSI-MI Controlled Vocabulary Terms

#### Interaction Detection Methods (MI:0001)

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0006 | anti bait coimmunoprecipitation | IP with antibody against bait |
| MI:0007 | anti tag coimmunoprecipitation | IP with antibody against tag |
| MI:0018 | two hybrid | Yeast two-hybrid |
| MI:0019 | coimmunoprecipitation | Co-IP |
| MI:0027 | cosedimentation | Co-sedimentation |
| MI:0045 | experimental interaction detection | Experimental method |
| MI:0071 | molecular sieving | Size exclusion |
| MI:0096 | pull down | Pull-down assay |
| MI:0114 | x-ray crystallography | Crystal structure |
| MI:0676 | tandem affinity purification | TAP tag purification |

#### Interaction Types (MI:0190)

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0195 | covalent binding | Covalent bond |
| MI:0407 | direct interaction | Direct physical contact |
| MI:0403 | colocalization | Same cellular location |
| MI:0914 | association | General association |
| MI:0915 | physical association | Physical contact |

#### Participant Biological Roles (MI:0500)

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0501 | enzyme | Catalytic activity |
| MI:0502 | enzyme target | Substrate |
| MI:0580 | electron donor | Reduces partner |
| MI:0581 | electron acceptor | Oxidizes partner |
| MI:0499 | unspecified role | Role not specified |

#### Participant Experimental Roles (MI:0495)

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0496 | bait | Tagged/immobilized protein |
| MI:0497 | neutral component | Neither bait nor prey |
| MI:0498 | prey | Co-purified protein |
| MI:0503 | self | Homodimer |

#### Interactor Types (MI:0313)

| MI ID | Term | Description |
|-------|------|-------------|
| MI:0326 | protein | Protein molecule |
| MI:0319 | dna | DNA molecule |
| MI:0320 | rna | RNA molecule |
| MI:0328 | small molecule | Chemical compound |
| MI:0324 | complex | Macromolecular complex |
| MI:2190 | gene | Gene locus |

---

## Confidence Scoring

### IntAct MI Score (intact-miscore)

IntAct provides a confidence score (0-1) based on:

| Factor | Weight | Description |
|--------|--------|-------------|
| Publication count | High | Multiple publications |
| Method diversity | High | Different experimental methods |
| Evidence type | Medium | Direct vs. indirect evidence |
| Curation depth | Medium | Level of annotation detail |

### Score Interpretation

| Score Range | Confidence | Description |
|-------------|------------|-------------|
| 0.75 - 1.00 | High | Multiple evidence types, well-curated |
| 0.50 - 0.74 | Medium | Good evidence, single publication |
| 0.25 - 0.49 | Low | Limited evidence |
| 0.00 - 0.24 | Very Low | Sparse annotation |

---

## Data Downloads

### FTP Structure

```
ftp://ftp.ebi.ac.uk/pub/databases/intact/current/
├── psimi/
│   ├── pmid/           # Per-publication files
│   ├── species/        # Per-species files
│   └── all/            # Complete dataset
├── psi25/              # PSI-MI XML 2.5 format
├── psi30/              # PSI-MI XML 3.0 format
└── various/            # Additional formats
```

### Key Download Files

| File | Description | Format |
|------|-------------|--------|
| `intact.txt` | Complete MITAB dataset | TSV |
| `intact.zip` | Complete PSI-MI XML | ZIP |
| `intact-micluster.txt` | Clustered interactions | TSV |
| `human.zip` | Human interactions only | XML |
| `mutation.txt` | Mutation annotations | TSV |

### Species-Specific Downloads

| Species | File | Size (approx) |
|---------|------|---------------|
| Homo sapiens | human.zip | 180 MB |
| Mus musculus | mouse.zip | 45 MB |
| Saccharomyces cerevisiae | yeast.zip | 120 MB |
| Drosophila melanogaster | fly.zip | 35 MB |

---

## Cross-References

### Supported Database Cross-References

| Database | ID Format | Example |
|----------|-----------|---------|
| UniProtKB | [A-Z][0-9]{5} or [A-Z][0-9][A-Z0-9]{3}[0-9] | P04637 |
| IntAct | EBI-####### | EBI-366083 |
| IMEx | IM-##### | IM-12345 |
| Ensembl | ENSG########### | ENSG00000141510 |
| RefSeq | NP_###### | NP_000537 |
| ChEBI | CHEBI:##### | CHEBI:15377 |
| Gene Ontology | GO:####### | GO:0005515 |
| PubMed | ######## | 1535557 |
| DOI | 10.####/... | 10.1016/0092-8674(92)90118-U |
| Reactome | R-HSA-###### | R-HSA-109581 |

### IMEx Consortium Partners

IntAct is part of the IMEx consortium with shared data from:

| Database | MI ID | Description |
|----------|-------|-------------|
| IntAct | MI:0469 | EMBL-EBI (primary) |
| MINT | MI:0471 | Molecular INTeraction database |
| DIP | MI:0465 | Database of Interacting Proteins |
| MatrixDB | MI:0917 | Extracellular matrix interactions |
| HPIDB | MI:0970 | Host-pathogen interactions |
| InnateDB | MI:0974 | Innate immunity interactions |
| MPIDB | MI:0903 | Microbial Protein Interactions |

---

## Sample Data

### Sample MITAB 2.5 Record

```
uniprotkb:P04637	uniprotkb:Q00987	intact:EBI-366083	intact:EBI-389768	uniprotkb:TP53(gene name)	uniprotkb:MDM2(gene name)	psi-mi:"MI:0018"(two hybrid)	Momand J (1992)	pubmed:1535557	taxid:9606(human)	taxid:9606(human)	psi-mi:"MI:0915"(physical association)	psi-mi:"MI:0469"(IntAct)	intact:EBI-77734	intact-miscore:0.73
```

### Sample JSON Response

```json
{
  "interactionAc": "EBI-77734",
  "interactorA": {
    "identifier": "P04637",
    "database": "uniprotkb",
    "shortLabel": "p53_human",
    "fullName": "Cellular tumor antigen p53",
    "aliases": [
      {"name": "TP53", "type": "gene name"}
    ],
    "organism": {
      "taxId": 9606,
      "scientificName": "Homo sapiens",
      "commonName": "Human"
    },
    "interactorType": {
      "miIdentifier": "MI:0326",
      "shortLabel": "protein"
    }
  },
  "interactorB": {
    "identifier": "Q00987",
    "database": "uniprotkb",
    "shortLabel": "mdm2_human",
    "fullName": "E3 ubiquitin-protein ligase Mdm2",
    "aliases": [
      {"name": "MDM2", "type": "gene name"}
    ],
    "organism": {
      "taxId": 9606,
      "scientificName": "Homo sapiens"
    },
    "interactorType": {
      "miIdentifier": "MI:0326",
      "shortLabel": "protein"
    }
  },
  "interactionType": {
    "miIdentifier": "MI:0915",
    "shortLabel": "physical association"
  },
  "detectionMethod": {
    "miIdentifier": "MI:0018",
    "shortLabel": "two hybrid"
  },
  "publication": {
    "pubmedId": "1535557",
    "firstAuthor": "Momand J",
    "year": 1992
  },
  "confidenceScore": 0.73,
  "sourceDatabase": "intact"
}
```

---

## Integration Notes

### Data Quality Indicators

| Indicator | Description |
|-----------|-------------|
| **IMEx curated** | Meets IMEx consortium standards |
| **Multiple publications** | Interaction reported in >1 paper |
| **Method diversity** | Detected by multiple methods |
| **Direct interaction** | MI:0407 type annotation |
| **Negative interaction** | Explicitly tested, no interaction |

### Best Practices for Integration

1. **Use UniProt identifiers** as primary keys for cross-referencing
2. **Filter by MI score** >= 0.5 for high-confidence interactions
3. **Consider detection method** - two-hybrid has higher false positive rate
4. **Check for complexes** - spoke vs matrix expansion affects interpretation
5. **Merge with IMEx partners** for comprehensive coverage

### API Considerations

- Rate limits apply; batch queries when possible
- Use MITAB format for large downloads
- PSI-MI XML provides richer annotation detail
- PSICQUIC allows federated queries across IMEx partners

---

## Data Format

| Format | Description | Use Case |
|--------|-------------|----------|
| Primary | MITAB 2.6 (Tab-separated) | Standard exchange format |
| Alternative | PSI-MI 3.0 XML | Rich annotation, standards-compliant |
| Alternative | JSON | Programmatic access |
| Search | PSICQUIC | Federated query protocol |
| Encoding | UTF-8 | All formats |

---

## Download

### Data Access Methods

| Method | URL | Format |
|--------|-----|--------|
| **Web Interface** | https://www.ebi.ac.uk/intact/home | Interactive search |
| **PSICQUIC Services** | https://www.ebi.ac.uk/Tools/webservices/psicquic | REST API |
| **FTP** | ftp://ftp.ebi.ac.uk/pub/databases/intact/ | MITAB, PSI-MI XML |
| **Downloads** | https://www.ebi.ac.uk/intact/download | Bulk exports |

### Download URLs

```bash
# MITAB26 format (full database)
ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip

# PSI-MI XML
ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psi30/intact.zip

# PSICQUIC query example
curl "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/query/query=*:*&format=tab25"
```

---

## Data Set Size

| Component | Records | Size (Est.) |
|-----------|---------|------------|
| **Protein Interactions** | 1,000,000+ | ~100-200 MB (MITAB) |
| **Binary Interactions** | 500,000+ (spoke) | ~50-100 MB |
| **Unique Proteins** | 500,000+ | Included in interactions |
| **Complexes** | 50,000+ | ~10 MB |
| **Total ZIP (MITAB)** | Full database | ~40-50 MB |
| **Total ZIP (XML)** | Full database | ~200-300 MB |

---

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `id` | string | Primary identifier | "EBI-1" |
| `name` | string | Entity name | "Protein A" |
| `type` | string | Record type | "protein" / "interaction" |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `interacts_with` | Protein | N:M |
| `described_in` | Publication | N:M |
| `has_participant` | Interactor | N:M |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `interactionAc` | IntAct interaction accession identifier | EBI-77734 |
| `intact-miscore` | IntAct molecular interaction confidence score (0-1) | 0.73 |
| `MI:XXXX` | PSI-MI controlled vocabulary term identifier | MI:0018 (two hybrid) |
| `taxid` | NCBI Taxonomy identifier for organism | taxid:9606 (human) |
| `bait` | Tagged/immobilized protein in pull-down experiments | MI:0496 |
| `prey` | Co-purified protein partner in interaction experiment | MI:0498 |
| `stoichiometry` | Number of molecules of each participant in complex | stoichiometry:2 |
| `rigid` | ROG-derived interaction digest checksum | rigid:4hL... |
| `rogid` | Relationship Ontology Graph ID for interactor checksum | rogid:UcdngwpTSS6hG... |
| `expansion method` | How complex data is converted to binary pairs | spoke, matrix |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| Binary Interaction | Pairwise association between two molecules | MITAB format |
| Physical Association | Direct or indirect physical contact between molecules | MI:0915 |
| Direct Interaction | Proven physical contact without intermediates | MI:0407 |
| Colocalization | Molecules found in same cellular location | MI:0403 |
| Detection Method | Experimental technique used to identify interaction | MI:0001 hierarchy |
| Participant Role | Biological function of molecule in interaction | enzyme, substrate |
| Experimental Role | Role assigned by experimental design | bait, prey, neutral |
| Interactor Type | Molecular class of participant | protein, DNA, RNA, small molecule |
| Complex Expansion | Converting n-ary complex to binary pairs | spoke (hub-based), matrix (all pairs) |
| Negative Interaction | Experimentally tested pair that does not interact | Column 36 in MITAB |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| IntAct | Interaction Database | EMBL-EBI molecular interaction database |
| PSI-MI | Proteomics Standards Initiative - Molecular Interaction | Standard for interaction data |
| MITAB | PSI-MI Tab-delimited | Tab-delimited exchange format |
| PSICQUIC | PSI Common Query Interface | Federated query web service |
| IMEx | International Molecular Exchange | Consortium for interaction curation |
| MIQL | Molecular Interaction Query Language | Query syntax for PSICQUIC |
| EBI | European Bioinformatics Institute | IntAct host institution |
| HUPO | Human Proteome Organization | PSI parent organization |
| TAP | Tandem Affinity Purification | Interaction detection method |
| Co-IP | Coimmunoprecipitation | Interaction detection method |
| MI | Molecular Interaction | PSI-MI ontology prefix |
| CC BY | Creative Commons Attribution | License type |

---

## References

1. Orchard S, et al. (2014) "The MIntAct project--IntAct as a common curation platform for 11 molecular interaction databases." Nucleic Acids Res. 42(Database issue):D358-63.

2. del-Toro N, et al. (2022) "The IntAct database: efficient access to fine-grained molecular interaction data." Nucleic Acids Res. 50(D1):D648-D657.

3. Kerrien S, et al. (2012) "The IntAct molecular interaction database in 2012." Nucleic Acids Res. 40(Database issue):D301-5.

4. HUPO-PSI MI Working Group: https://github.com/HUPO-PSI/miXML

5. IMEx Consortium: https://www.imexconsortium.org/

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial schema documentation for IntAct |
