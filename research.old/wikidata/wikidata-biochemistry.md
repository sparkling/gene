# Wikidata Biochemistry Coverage

## Overview

Wikidata serves as a comprehensive semantic framework for biochemical data, integrating information from major biological databases including NCBI, UniProt, Reactome, KEGG, WikiPathways, ChEBI, and Gene Ontology. This document details Wikidata's coverage of biochemical pathways, genes, proteins, enzymes, and chemical compounds, along with comprehensive SPARQL queries for data extraction.

**Key Statistics:**
- ~59,721 human genes imported from NCBI
- ~27,306 human proteins imported from UniProt (SwissProt subset)
- ~3,000 human biological pathways (from Reactome and WikiPathways)
- ~150,000+ chemical compounds including drugs

**References:**
- [Wikidata as a semantic framework for the Gene Wiki initiative](https://ncbi.nlm.nih.gov/pmc/articles/PMC4795929)
- [Wikidata as a knowledge graph for the life sciences](https://pmc.ncbi.nlm.nih.gov/articles/PMC7077981/)

---

## 1. Biochemical Pathways

### Core Wikidata Items

| Item | Q-ID | Description |
|------|------|-------------|
| biological pathway | Q4915012 | Generic concept for biological pathways |
| biological process | Q2996394 | GO term for biological processes |
| metabolic pathway | Q4915059 | Pathways involving metabolism |
| signal transduction | Q189093 | Signaling pathways |

### Pathway Database Properties

| Property | ID | Description | Example Values |
|----------|-----|-------------|----------------|
| Reactome ID | P3937 | Reactome pathway identifier | R-HSA-71403 (TCA cycle) |
| KEGG ID | P665 | KEGG pathway/compound ID | hsa00010 (Glycolysis) |
| WikiPathways ID | P2410 | WikiPathways identifier | WP78 (TCA Cycle) |

### Pathway Data Model

Wikidata pathways contain:
- **Name and description** (multilingual labels)
- **Organism** (P703 - found in taxon)
- **External identifiers** (Reactome, KEGG, WikiPathways)
- **Component genes** (P527 - has part)
- **Component proteins** (P527)
- **Component compounds** (P527)

### SPARQL Queries for Pathways

#### Query 1: All Human Pathways with External Identifiers

```sparql
# All human biological pathways with Reactome, KEGG, and WikiPathways IDs
SELECT DISTINCT ?pathway ?pathwayLabel ?reactomeId ?keggId ?wikipathwaysId
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .  # instance of biological pathway (or subclass)
  ?pathway wdt:P703 wd:Q15978631 .           # found in taxon: Homo sapiens

  OPTIONAL { ?pathway wdt:P3937 ?reactomeId . }      # Reactome ID
  OPTIONAL { ?pathway wdt:P665 ?keggId . }           # KEGG ID
  OPTIONAL { ?pathway wdt:P2410 ?wikipathwaysId . }  # WikiPathways ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
```

#### Query 2: Metabolic Pathways with Component Genes

```sparql
# Human metabolic pathways with their component genes
SELECT DISTINCT ?pathway ?pathwayLabel ?gene ?geneLabel ?entrezId
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q4915059 .  # instance of metabolic pathway
  ?pathway wdt:P703 wd:Q15978631 .           # found in Homo sapiens
  ?pathway wdt:P527 ?gene .                  # has part (gene)
  ?gene wdt:P31 wd:Q7187 .                   # gene is instance of gene
  ?gene wdt:P351 ?entrezId .                 # Entrez Gene ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel ?geneLabel
LIMIT 1000
```

#### Query 3: Signal Transduction Pathways

```sparql
# Signal transduction pathways with associated proteins
SELECT DISTINCT ?pathway ?pathwayLabel ?protein ?proteinLabel ?uniprotId
WHERE {
  ?pathway wdt:P31/wdt:P279* wd:Q189093 .   # instance of signal transduction
  ?pathway wdt:P703 wd:Q15978631 .           # found in Homo sapiens
  ?pathway wdt:P527 ?protein .               # has part
  ?protein wdt:P31 wd:Q8054 .                # protein
  ?protein wdt:P352 ?uniprotId .             # UniProt ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
LIMIT 500
```

#### Query 4: Pathways by Reactome ID

```sparql
# Get pathway details by Reactome identifier
SELECT DISTINCT ?pathway ?pathwayLabel ?description ?organism ?organismLabel
WHERE {
  ?pathway wdt:P3937 ?reactomeId .           # has Reactome ID
  ?pathway wdt:P703 ?organism .              # found in taxon
  OPTIONAL { ?pathway schema:description ?description . FILTER(LANG(?description) = "en") }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
```

---

## 2. Genes

### Core Wikidata Items

| Item | Q-ID | Description |
|------|------|-------------|
| gene | Q7187 | Generic gene concept |
| protein-coding gene | Q20747295 | Gene that encodes protein |
| ncRNA gene | Q27087 | Non-coding RNA gene |
| pseudogene | Q277338 | Non-functional gene copy |

### Gene Properties

| Property | ID | Description |
|----------|-----|-------------|
| Entrez Gene ID | P351 | NCBI Gene identifier |
| HGNC gene symbol | P353 | Official gene symbol |
| HGNC ID | P354 | HGNC identifier |
| Ensembl gene ID | P594 | Ensembl gene identifier |
| Ensembl transcript ID | P704 | Ensembl transcript identifier |
| chromosome | P1057 | Chromosome location |
| genomic start | P644 | Start position |
| genomic end | P645 | End position |
| encodes | P688 | Protein encoded by gene |
| genetic association | P2293 | Gene-disease association |

### Human Gene Coverage

Wikidata contains ~59,721 human genes imported from NCBI, compared to:
- ~20,000 protein-coding genes
- ~25,000 non-coding genes
- Additional pseudogenes and other gene types

Gene data includes:
- Genomic coordinates (GRCh37 and GRCh38 assemblies)
- Cross-references to Entrez, HGNC, Ensembl
- Links to encoded proteins
- Disease associations

### SPARQL Queries for Genes

#### Query 5: All Human Genes with Identifiers

```sparql
# All human protein-coding genes with major identifiers
SELECT DISTINCT ?gene ?geneLabel ?entrezId ?hgncSymbol ?hgncId ?ensemblId
WHERE {
  ?gene wdt:P31 wd:Q20747295 .              # instance of protein-coding gene
  ?gene wdt:P703 wd:Q15978631 .             # found in taxon: Homo sapiens

  OPTIONAL { ?gene wdt:P351 ?entrezId . }    # Entrez Gene ID
  OPTIONAL { ?gene wdt:P353 ?hgncSymbol . }  # HGNC gene symbol
  OPTIONAL { ?gene wdt:P354 ?hgncId . }      # HGNC ID
  OPTIONAL { ?gene wdt:P594 ?ensemblId . }   # Ensembl gene ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?hgncSymbol
```

#### Query 6: Genes with Chromosomal Location

```sparql
# Human genes with chromosome and position (GRCh38)
SELECT DISTINCT ?gene ?geneLabel ?hgncSymbol ?chromosome ?chromosomeLabel ?start ?end
WHERE {
  ?gene wdt:P31/wdt:P279* wd:Q7187 .        # instance of gene
  ?gene wdt:P703 wd:Q15978631 .              # Homo sapiens
  ?gene wdt:P353 ?hgncSymbol .               # HGNC symbol

  ?gene p:P1057 ?chromStatement .            # chromosome statement
  ?chromStatement ps:P1057 ?chromosome .     # chromosome value
  ?chromStatement pq:P659 wd:Q20966585 .     # qualifier: GRCh38

  OPTIONAL {
    ?gene p:P644 ?startStatement .
    ?startStatement ps:P644 ?start .
    ?startStatement pq:P659 wd:Q20966585 .
  }
  OPTIONAL {
    ?gene p:P645 ?endStatement .
    ?endStatement ps:P645 ?end .
    ?endStatement pq:P659 wd:Q20966585 .
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?chromosome ?start
LIMIT 1000
```

#### Query 7: Gene-Disease Associations

```sparql
# Genes associated with diseases (using P2293 genetic association)
SELECT DISTINCT ?gene ?geneLabel ?hgncSymbol ?disease ?diseaseLabel ?diseaseOmimId
WHERE {
  ?gene wdt:P31/wdt:P279* wd:Q7187 .        # instance of gene
  ?gene wdt:P703 wd:Q15978631 .              # Homo sapiens
  ?gene wdt:P353 ?hgncSymbol .               # HGNC symbol

  ?disease wdt:P2293 ?gene .                 # genetic association to gene
  OPTIONAL { ?disease wdt:P492 ?diseaseOmimId . }  # OMIM ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?diseaseLabel ?geneLabel
LIMIT 1000
```

#### Query 8: Genes and Their Encoded Proteins

```sparql
# Human genes with their encoded proteins and UniProt IDs
SELECT DISTINCT ?gene ?geneLabel ?hgncSymbol ?protein ?proteinLabel ?uniprotId
WHERE {
  ?gene wdt:P31/wdt:P279* wd:Q7187 .        # gene
  ?gene wdt:P703 wd:Q15978631 .              # Homo sapiens
  ?gene wdt:P353 ?hgncSymbol .               # HGNC symbol
  ?gene wdt:P688 ?protein .                  # encodes protein
  ?protein wdt:P352 ?uniprotId .             # UniProt ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?hgncSymbol
```

---

## 3. Proteins

### Core Wikidata Items

| Item | Q-ID | Description |
|------|------|-------------|
| protein | Q8054 | Generic protein concept |
| enzyme | Q8047 | Catalytic proteins |
| receptor | Q193053 | Cell receptors |
| transcription factor | Q3271540 | DNA-binding proteins |

### Protein Properties

| Property | ID | Description |
|----------|-----|-------------|
| UniProt protein ID | P352 | UniProt identifier |
| Ensembl protein ID | P705 | Ensembl protein identifier |
| encoded by | P702 | Gene that encodes protein |
| molecular function | P680 | GO molecular function |
| cell component | P681 | GO cellular component |
| biological process | P682 | GO biological process |
| EC number | P591 | Enzyme Commission number |
| found in taxon | P703 | Organism |

### Human Protein Coverage

Wikidata contains ~27,306 human proteins imported from the SwissProt subset of UniProt, with:
- UniProt identifiers
- Ensembl protein IDs
- Gene Ontology annotations (molecular function, cellular component, biological process)
- Links to encoding genes
- Enzyme classifications (for enzymes)

### SPARQL Queries for Proteins

#### Query 9: All Human Proteins with Identifiers

```sparql
# All human proteins with UniProt and Ensembl IDs
SELECT DISTINCT ?protein ?proteinLabel ?uniprotId ?ensemblProteinId ?gene ?geneLabel
WHERE {
  ?protein wdt:P31 wd:Q8054 .               # instance of protein
  ?protein wdt:P703 wd:Q15978631 .          # found in Homo sapiens
  ?protein wdt:P352 ?uniprotId .            # UniProt ID (required)

  OPTIONAL { ?protein wdt:P705 ?ensemblProteinId . }  # Ensembl protein ID
  OPTIONAL { ?protein wdt:P702 ?gene . }              # encoded by gene

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?proteinLabel
```

#### Query 10: Proteins with GO Annotations

```sparql
# Human proteins with Gene Ontology annotations
SELECT DISTINCT ?protein ?proteinLabel ?uniprotId
       ?molecularFunction ?molecularFunctionLabel
       ?cellComponent ?cellComponentLabel
       ?biologicalProcess ?biologicalProcessLabel
WHERE {
  ?protein wdt:P31 wd:Q8054 .               # protein
  ?protein wdt:P703 wd:Q15978631 .          # Homo sapiens
  ?protein wdt:P352 ?uniprotId .            # UniProt ID

  OPTIONAL { ?protein wdt:P680 ?molecularFunction . }   # GO molecular function
  OPTIONAL { ?protein wdt:P681 ?cellComponent . }       # GO cellular component
  OPTIONAL { ?protein wdt:P682 ?biologicalProcess . }   # GO biological process

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?proteinLabel
LIMIT 500
```

#### Query 11: Proteins by GO Term

```sparql
# Find all human proteins with a specific GO molecular function
# Example: kinase activity (GO:0016301 = Q21099621)
SELECT DISTINCT ?protein ?proteinLabel ?uniprotId ?goTerm ?goTermLabel
WHERE {
  VALUES ?goTerm { wd:Q21099621 }           # kinase activity

  ?protein wdt:P31 wd:Q8054 .               # protein
  ?protein wdt:P703 wd:Q15978631 .          # Homo sapiens
  ?protein wdt:P352 ?uniprotId .            # UniProt ID
  ?protein wdt:P680 ?goTerm .               # has GO molecular function

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?proteinLabel
```

#### Query 12: Membrane Proteins (Cellular Component)

```sparql
# Human proteins located in plasma membrane
SELECT DISTINCT ?protein ?proteinLabel ?uniprotId
WHERE {
  ?protein wdt:P31 wd:Q8054 .               # protein
  ?protein wdt:P703 wd:Q15978631 .          # Homo sapiens
  ?protein wdt:P352 ?uniprotId .            # UniProt ID
  ?protein wdt:P681 wd:Q74615596 .          # cellular component: plasma membrane

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?proteinLabel
```

---

## 4. Enzymes

### Core Wikidata Items

| Item | Q-ID | Description |
|------|------|-------------|
| enzyme | Q8047 | Catalytic protein |
| oxidoreductase | Q178211 | EC 1.x.x.x |
| transferase | Q178401 | EC 2.x.x.x |
| hydrolase | Q182858 | EC 3.x.x.x |
| lyase | Q180891 | EC 4.x.x.x |
| isomerase | Q182847 | EC 5.x.x.x |
| ligase | Q178352 | EC 6.x.x.x |
| cytochrome P450 | Q423111 | CYP450 enzyme superfamily |

### Enzyme Properties

| Property | ID | Description |
|----------|-----|-------------|
| EC number | P591 | Enzyme Commission classification |
| found in taxon | P703 | Organism source |
| catalyzes | P3923 | Reaction catalyzed (links to Rhea) |
| cofactor | P5974 | Required cofactors |

### Pharmacogenomically Relevant Enzymes

Cytochrome P450 (CYP) enzymes are particularly important for drug metabolism:
- CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP3A4 (major drug-metabolizing enzymes)
- CYP nomenclature: Family number + subfamily letter + individual enzyme number (e.g., CYP2D6)

### SPARQL Queries for Enzymes

#### Query 13: All Human Enzymes with EC Numbers

```sparql
# Human enzymes with EC classification
SELECT DISTINCT ?enzyme ?enzymeLabel ?ecNumber ?uniprotId
WHERE {
  ?enzyme wdt:P31/wdt:P279* wd:Q8047 .      # instance of enzyme (or subclass)
  ?enzyme wdt:P703 wd:Q15978631 .            # Homo sapiens
  ?enzyme wdt:P591 ?ecNumber .               # EC number

  OPTIONAL { ?enzyme wdt:P352 ?uniprotId . } # UniProt ID

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?ecNumber
```

#### Query 14: Cytochrome P450 Enzymes

```sparql
# CYP450 enzymes in humans (pharmacogenomics relevant)
SELECT DISTINCT ?enzyme ?enzymeLabel ?ecNumber ?uniprotId ?geneLabel ?hgncSymbol
WHERE {
  ?enzyme wdt:P31/wdt:P279* wd:Q423111 .    # subclass of cytochrome P450
  ?enzyme wdt:P703 wd:Q15978631 .            # Homo sapiens

  OPTIONAL { ?enzyme wdt:P591 ?ecNumber . }  # EC number
  OPTIONAL { ?enzyme wdt:P352 ?uniprotId . } # UniProt ID
  OPTIONAL {
    ?enzyme wdt:P702 ?gene .                 # encoded by
    ?gene wdt:P353 ?hgncSymbol .             # HGNC symbol
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?enzymeLabel
```

#### Query 15: Enzymes by EC Class

```sparql
# Find all human oxidoreductases (EC 1.x.x.x)
SELECT DISTINCT ?enzyme ?enzymeLabel ?ecNumber ?uniprotId
WHERE {
  ?enzyme wdt:P31/wdt:P279* wd:Q178211 .    # oxidoreductase
  ?enzyme wdt:P703 wd:Q15978631 .            # Homo sapiens

  OPTIONAL { ?enzyme wdt:P591 ?ecNumber . }
  OPTIONAL { ?enzyme wdt:P352 ?uniprotId . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?ecNumber
LIMIT 500
```

#### Query 16: Enzyme-Substrate Relationships

```sparql
# Enzymes and their substrates (via physically interacts with + role qualifier)
SELECT DISTINCT ?enzyme ?enzymeLabel ?substrate ?substrateLabel ?role ?roleLabel
WHERE {
  ?enzyme wdt:P31/wdt:P279* wd:Q8047 .      # enzyme
  ?enzyme wdt:P703 wd:Q15978631 .            # Homo sapiens

  ?enzyme p:P129 ?interaction .              # physically interacts with
  ?interaction ps:P129 ?substrate .          # substrate
  ?interaction pq:P2868 ?role .              # subject has role

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?enzymeLabel
LIMIT 500
```

---

## 5. Chemical Compounds (Biochemical)

### Core Wikidata Items

| Item | Q-ID | Description |
|------|------|-------------|
| chemical compound | Q11173 | Generic chemical compound |
| metabolite | Q25617064 | Metabolic intermediate |
| cofactor | Q28631133 | Enzyme cofactor |
| hormone | Q11364 | Signaling molecule |
| vitamin | Q34956 | Essential micronutrient |
| amino acid | Q8066 | Protein building blocks |
| nucleotide | Q27242 | DNA/RNA building blocks |

### Chemical Compound Properties

| Property | ID | Description |
|----------|-----|-------------|
| ChEBI ID | P683 | ChEBI identifier |
| PubChem CID | P662 | PubChem Compound ID |
| CAS Registry Number | P231 | CAS identifier |
| ECHA InfoCard ID | P2566 | European chemical identifier |
| InChIKey | P235 | Chemical hash |
| SMILES | P233 | Structure notation |
| physically interacts with | P129 | Compound-protein interactions |
| KEGG ID | P665 | KEGG compound ID |

### Compound Data Sources

Wikidata integrates chemical data from:
- ChEBI (Chemical Entities of Biological Interest) - ~195,000 entries
- PubChem - millions of compounds
- KEGG Compound database
- HMDB (Human Metabolome Database)
- DrugBank

### SPARQL Queries for Chemical Compounds

#### Query 17: Metabolites with Database IDs

```sparql
# Metabolites with ChEBI, PubChem, and KEGG identifiers
SELECT DISTINCT ?compound ?compoundLabel ?chebiId ?pubchemCid ?keggId ?casNumber
WHERE {
  ?compound wdt:P31/wdt:P279* wd:Q25617064 . # metabolite

  OPTIONAL { ?compound wdt:P683 ?chebiId . }     # ChEBI ID
  OPTIONAL { ?compound wdt:P662 ?pubchemCid . }  # PubChem CID
  OPTIONAL { ?compound wdt:P665 ?keggId . }      # KEGG ID
  OPTIONAL { ?compound wdt:P231 ?casNumber . }   # CAS number

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?compoundLabel
LIMIT 1000
```

#### Query 18: Hormones and Their Targets

```sparql
# Human hormones and their protein targets
SELECT DISTINCT ?hormone ?hormoneLabel ?target ?targetLabel ?uniprotId
WHERE {
  ?hormone wdt:P31/wdt:P279* wd:Q11364 .    # hormone
  ?hormone wdt:P129 ?target .                # physically interacts with
  ?target wdt:P31 wd:Q8054 .                 # target is protein
  ?target wdt:P703 wd:Q15978631 .            # human protein

  OPTIONAL { ?target wdt:P352 ?uniprotId . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?hormoneLabel
```

#### Query 19: Enzyme Cofactors

```sparql
# Cofactors and their associated enzymes
SELECT DISTINCT ?cofactor ?cofactorLabel ?chebiId ?enzyme ?enzymeLabel
WHERE {
  ?cofactor wdt:P31/wdt:P279* wd:Q28631133 . # cofactor

  OPTIONAL { ?cofactor wdt:P683 ?chebiId . }

  # Find enzymes that use this cofactor
  OPTIONAL {
    ?enzyme wdt:P5974 ?cofactor .             # has cofactor
    ?enzyme wdt:P703 wd:Q15978631 .           # human enzyme
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?cofactorLabel
LIMIT 500
```

#### Query 20: Compounds by ChEBI ID

```sparql
# Get compound details by ChEBI identifier
SELECT DISTINCT ?compound ?compoundLabel ?formula ?inchikey ?smiles
WHERE {
  ?compound wdt:P683 "CHEBI:15377" .         # water as example

  OPTIONAL { ?compound wdt:P274 ?formula . }   # molecular formula
  OPTIONAL { ?compound wdt:P235 ?inchikey . }  # InChIKey
  OPTIONAL { ?compound wdt:P233 ?smiles . }    # canonical SMILES

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
```

---

## 6. Drug-Gene-Pathway Links

### Traversal Pattern: Drug -> Target Protein -> Pathway

The key traversal properties:
1. **Drug -> Protein**: `wdt:P129` (physically interacts with)
2. **Protein -> Gene**: `wdt:P702` (encoded by)
3. **Gene -> Disease**: `wdt:P2293` (genetic association)
4. **Pathway -> Components**: `wdt:P527` (has part)

### SPARQL Queries for Integrated Analysis

#### Query 21: Drug Targets and Associated Pathways

```sparql
# Drugs, their protein targets, and pathways containing those targets
SELECT DISTINCT ?drug ?drugLabel ?target ?targetLabel ?pathway ?pathwayLabel
WHERE {
  ?drug wdt:P31 wd:Q12140 .                  # instance of medication
  ?drug wdt:P129 ?target .                   # drug interacts with target
  ?target wdt:P31 wd:Q8054 .                 # target is protein
  ?target wdt:P703 wd:Q15978631 .            # human protein

  # Find pathways containing this target
  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .   # biological pathway
  ?pathway wdt:P527 ?target .                 # pathway has part target

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?drugLabel
LIMIT 500
```

#### Query 22: Drug -> Target -> Gene -> Disease

```sparql
# Complete drug-target-gene-disease chain
SELECT DISTINCT ?drug ?drugLabel ?protein ?proteinLabel
       ?gene ?geneLabel ?disease ?diseaseLabel
WHERE {
  ?drug wdt:P31 wd:Q12140 .                  # medication
  ?drug wdt:P129 ?protein .                  # interacts with protein
  ?protein wdt:P31 wd:Q8054 .                # protein
  ?protein wdt:P703 wd:Q15978631 .           # human
  ?protein wdt:P702 ?gene .                  # encoded by gene
  ?disease wdt:P2293 ?gene .                 # disease associated with gene

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?drugLabel
LIMIT 500
```

#### Query 23: Cancer Drugs and Biological Processes

```sparql
# Drugs targeting cancer pathways with biological process annotations
SELECT DISTINCT ?drug ?drugLabel ?protein ?proteinLabel
       ?bioProcess ?bioProcessLabel ?disease ?diseaseLabel
WHERE {
  VALUES ?cancerGoTerms { wd:Q14818032 wd:Q14599311 wd:Q2383867 }  # cell proliferation terms

  ?drug wdt:P129 ?protein .                  # drug interacts with protein
  ?gene wdt:P688 ?protein .                  # gene encodes protein
  ?disease wdt:P2293 ?gene .                 # genetic association
  ?disease wdt:P279* wd:Q12078 .             # subclass of cancer
  ?protein wdt:P682 ?bioProcess .            # biological process
  ?bioProcess (wdt:P361|wdt:P279)* ?cancerGoTerms .  # related to cancer processes

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 200
```

---

## 7. Federated SPARQL Queries

### Available SPARQL Endpoints for Federation

| Endpoint | URL | Description |
|----------|-----|-------------|
| Wikidata | https://query.wikidata.org/sparql | Main Wikidata endpoint |
| UniProt | https://sparql.uniprot.org/sparql | Protein data |
| WikiPathways | https://sparql.wikipathways.org/sparql | Pathway data |
| ChEBI | https://www.ebi.ac.uk/rdf/services/sparql | Chemical entities |
| Rhea | https://sparql.rhea-db.org/sparql | Biochemical reactions |
| neXtProt | https://sparql.nextprot.org/sparql | Human protein data |

### Federated Query Examples

#### Query 24: Wikidata + UniProt Federation

```sparql
# Federated query: Wikidata drugs that interact with UniProt proteins
# Run this on the UniProt SPARQL endpoint
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wd: <http://www.wikidata.org/entity/>

SELECT DISTINCT ?drug ?drugLabel ?uniprotEntry ?proteinName
WHERE {
  # Wikidata part
  SERVICE <https://query.wikidata.org/sparql> {
    ?drug wdt:P31 wd:Q12140 .               # medication
    ?drug wdt:P129 ?target .                # interacts with
    ?target wdt:P352 ?uniprotId .           # UniProt ID
    ?drug rdfs:label ?drugLabel .
    FILTER(LANG(?drugLabel) = "en")
  }

  # UniProt part
  BIND(IRI(CONCAT("http://purl.uniprot.org/uniprot/", ?uniprotId)) AS ?uniprotEntry)
  ?uniprotEntry up:recommendedName ?recName .
  ?recName up:fullName ?proteinName .
}
LIMIT 100
```

#### Query 25: Wikidata + WikiPathways Federation

```sparql
# Federated query combining Wikidata metabolites with WikiPathways
# Run on WikiPathways endpoint
PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX wdt: <http://www.wikidata.org/prop/direct/>
PREFIX wd: <http://www.wikidata.org/entity/>

SELECT DISTINCT ?metabolite ?metaboliteLabel ?pathway ?pathwayTitle
WHERE {
  # WikiPathways part
  ?pathway a wp:Pathway .
  ?pathway wp:organismName "Homo sapiens" .
  ?pathway dc:title ?pathwayTitle .
  ?pathway wp:isAbout ?pwDataNode .
  ?pwDataNode wp:bdbWikidata ?metabolite .

  # Wikidata part for labels
  SERVICE <https://query.wikidata.org/sparql> {
    ?metabolite rdfs:label ?metaboliteLabel .
    FILTER(LANG(?metaboliteLabel) = "en")
  }
}
LIMIT 200
```

#### Query 26: Enzymes and Rhea Reactions

```sparql
# Find Rhea reactions catalyzed by human enzymes
# Run on Wikidata endpoint
PREFIX rh: <http://rdf.rhea-db.org/>

SELECT DISTINCT ?enzyme ?enzymeLabel ?rheaReaction ?ecNumber
WHERE {
  ?enzyme wdt:P31/wdt:P279* wd:Q8047 .      # enzyme
  ?enzyme wdt:P703 wd:Q15978631 .            # human
  ?enzyme wdt:P591 ?ecNumber .               # EC number

  # Get Rhea reaction via federated query
  SERVICE <https://sparql.rhea-db.org/sparql> {
    ?rheaReaction rh:ec ?ecNumberUri .
    FILTER(CONTAINS(STR(?ecNumberUri), ?ecNumber))
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 100
```

---

## 8. Practical Use Cases

### Use Case 1: Pharmacogenomics Pipeline

```sparql
# Find drugs that interact with pharmacogenomically important CYP enzymes
SELECT DISTINCT ?drug ?drugLabel ?enzyme ?enzymeLabel ?gene ?geneSymbol
WHERE {
  VALUES ?cyp { wd:Q423111 }                 # CYP450 superfamily

  ?drug wdt:P31 wd:Q12140 .                  # medication
  ?drug wdt:P129 ?enzyme .                   # interacts with enzyme
  ?enzyme wdt:P31/wdt:P279* ?cyp .           # enzyme is CYP450
  ?enzyme wdt:P703 wd:Q15978631 .            # human

  OPTIONAL {
    ?enzyme wdt:P702 ?gene .
    ?gene wdt:P353 ?geneSymbol .
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?drugLabel
```

### Use Case 2: Disease-Pathway Analysis

```sparql
# Find pathways relevant to a specific disease
# Example: Diabetes mellitus (Q12206)
SELECT DISTINCT ?pathway ?pathwayLabel ?gene ?geneLabel ?protein ?proteinLabel
WHERE {
  wd:Q12206 wdt:P2293 ?gene .                # genes associated with diabetes
  ?gene wdt:P688 ?protein .                  # gene encodes protein

  ?pathway wdt:P31/wdt:P279* wd:Q4915012 .   # biological pathway
  ?pathway wdt:P527 ?protein .               # pathway contains protein

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?pathwayLabel
```

### Use Case 3: Metabolite-Enzyme Network

```sparql
# Build network of metabolites and enzymes that process them
SELECT DISTINCT ?metabolite ?metaboliteLabel ?enzyme ?enzymeLabel ?ecNumber
WHERE {
  ?metabolite wdt:P31/wdt:P279* wd:Q25617064 .  # metabolite
  ?metabolite wdt:P129 ?enzyme .                 # interacts with enzyme
  ?enzyme wdt:P31/wdt:P279* wd:Q8047 .           # enzyme
  ?enzyme wdt:P703 wd:Q15978631 .                # human

  OPTIONAL { ?enzyme wdt:P591 ?ecNumber . }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
ORDER BY ?metaboliteLabel
LIMIT 500
```

---

## 9. Data Quality and Maintenance

### Consistency Checks

Wikidata uses SPARQL-based consistency tests to detect data model issues:

```sparql
# Check for genes without species annotation
SELECT ?gene ?geneLabel WHERE {
  ?gene wdt:P31 wd:Q7187 .                   # gene
  FILTER NOT EXISTS { ?gene wdt:P703 ?taxon . }  # no taxon
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 100
```

```sparql
# Check for proteins without encoding gene
SELECT ?protein ?proteinLabel WHERE {
  ?protein wdt:P31 wd:Q8054 .                # protein
  ?protein wdt:P703 wd:Q15978631 .           # human
  FILTER NOT EXISTS { ?protein wdt:P702 ?gene . }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE],en". }
}
LIMIT 100
```

### Bot-Maintained Data

Wikidata biological data is maintained by automated bots including:
- **ProteinBoxBot**: Imports gene and protein data from NCBI and UniProt
- **GeneWiki bots**: Maintain gene-disease associations
- **WikiPathways bot**: Syncs pathway data

---

## 10. Summary: Key Property Reference

### Essential Biological Properties

| Category | Property | ID | Usage |
|----------|----------|-----|-------|
| **Classification** | instance of | P31 | Item type classification |
| | subclass of | P279 | Hierarchical relationships |
| **Organism** | found in taxon | P703 | Species association |
| **Gene IDs** | Entrez Gene ID | P351 | NCBI gene identifier |
| | HGNC gene symbol | P353 | Official gene symbol |
| | Ensembl gene ID | P594 | Ensembl identifier |
| **Protein IDs** | UniProt protein ID | P352 | UniProt identifier |
| | Ensembl protein ID | P705 | Ensembl protein ID |
| **Gene Ontology** | molecular function | P680 | GO MF annotation |
| | cell component | P681 | GO CC annotation |
| | biological process | P682 | GO BP annotation |
| **Pathways** | Reactome ID | P3937 | Reactome identifier |
| | KEGG ID | P665 | KEGG identifier |
| | WikiPathways ID | P2410 | WikiPathways ID |
| **Chemistry** | ChEBI ID | P683 | ChEBI identifier |
| | PubChem CID | P662 | PubChem compound ID |
| | CAS number | P231 | CAS registry number |
| **Relationships** | encodes | P688 | Gene to protein |
| | encoded by | P702 | Protein to gene |
| | genetic association | P2293 | Gene-disease link |
| | physically interacts | P129 | Molecular interactions |
| | EC number | P591 | Enzyme classification |

---

## References

1. [Wikidata Query Service](https://query.wikidata.org/) - Main SPARQL endpoint
2. [Wikidata SPARQL Tutorial](https://www.wikidata.org/wiki/Wikidata:SPARQL_tutorial)
3. [Wikidata SPARQL Examples](https://www.wikidata.org/wiki/Wikidata:SPARQL_query_service/queries/examples)
4. [Gene Wiki Initiative](https://ncbi.nlm.nih.gov/pmc/articles/PMC4795929)
5. [WikiPathways SPARQL](https://www.wikipathways.org/sparql.html)
6. [UniProt SPARQL](https://sparql.uniprot.org/)
7. [ChEBI Database](https://www.ebi.ac.uk/chebi/)
8. [Reactome Pathway Database](https://reactome.org/)
9. [KEGG Pathway Database](https://www.genome.jp/kegg/pathway.html)
10. [HGNC Gene Nomenclature](https://www.genenames.org/)

---

*Document generated: 2025*
*Wikidata SPARQL endpoint: https://query.wikidata.org/sparql*
