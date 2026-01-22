# Pathway Database Sample Data

This document contains actual sample data retrieved from Reactome, WikiPathways, and DisGeNET APIs.

---

## 1. Reactome Sample Data

### Pathway: Metabolism (R-HSA-1430728)

**API Call:** `GET https://reactome.org/ContentService/data/query/R-HSA-1430728`

```json
{
  "dbId": 1430728,
  "displayName": "Metabolism",
  "stId": "R-HSA-1430728",
  "stIdVersion": "R-HSA-1430728.16",
  "isInDisease": false,
  "isInferred": false,
  "maxDepth": 7,
  "name": ["Metabolism"],
  "releaseDate": "2011-09-20",
  "releaseStatus": "UPDATED",
  "speciesName": "Homo sapiens",
  "goBiologicalProcess": {
    "dbId": 504,
    "displayName": "metabolic process",
    "accession": "0008152",
    "databaseName": "GO",
    "url": "https://www.ebi.ac.uk/QuickGO/term/GO:0008152"
  },
  "species": [{
    "dbId": 48887,
    "displayName": "Homo sapiens",
    "taxId": "9606",
    "abbreviation": "HSA"
  }],
  "doi": "10.3180/R-HSA-1430728.15",
  "hasDiagram": true,
  "hasEHLD": true,
  "schemaClass": "TopLevelPathway"
}
```

### Pathway: Cholesterol Biosynthesis (R-HSA-191273)

**API Call:** `GET https://reactome.org/ContentService/data/query/R-HSA-191273`

```json
{
  "dbId": 191273,
  "displayName": "Cholesterol biosynthesis",
  "stId": "R-HSA-191273",
  "stIdVersion": "R-HSA-191273.13",
  "isInDisease": false,
  "maxDepth": 3,
  "speciesName": "Homo sapiens",
  "compartment": [
    {
      "dbId": 70101,
      "displayName": "cytosol",
      "accession": "0005829",
      "databaseName": "GO"
    },
    {
      "dbId": 12045,
      "displayName": "endoplasmic reticulum membrane",
      "accession": "0005789",
      "databaseName": "GO"
    }
  ],
  "goBiologicalProcess": {
    "dbId": 20521,
    "displayName": "cholesterol biosynthetic process",
    "accession": "0006695",
    "databaseName": "GO"
  },
  "literatureReference": [
    {
      "dbId": 194670,
      "title": "Membrane-bound enzymes of cholesterol synthesis from lanosterol",
      "journal": "Biochem Biophys Res Commun",
      "pubMedIdentifier": 11969204,
      "year": 2002
    }
  ],
  "doi": "10.3180/R-HSA-191273.11",
  "hasDiagram": true,
  "hasEHLD": false,
  "schemaClass": "Pathway"
}
```

### Small Molecule: ATP (R-ALL-29358)

**API Call:** `GET https://reactome.org/ContentService/data/query/29358`

```json
{
  "dbId": 29358,
  "displayName": "ATP [nucleoplasm]",
  "stId": "R-ALL-29358",
  "stIdVersion": "R-ALL-29358.3",
  "maxDepth": 1,
  "name": ["ATP", "Adenosine 5'-triphosphate", "ATP(4-)"],
  "compartment": [{
    "dbId": 7660,
    "displayName": "nucleoplasm",
    "accession": "0005654",
    "databaseName": "GO"
  }],
  "referenceEntity": {
    "dbId": 8869364,
    "stId": "chebi:30616",
    "databaseName": "ChEBI",
    "identifier": "30616",
    "name": ["ATP(4-)", "ATP", "atp", "Adenosine 5'-triphosphate"],
    "formula": "C10H12N5O13P3",
    "schemaClass": "ReferenceMolecule"
  },
  "schemaClass": "SimpleEntity"
}
```

### Cholesterol Biosynthesis Participants

**API Call:** `GET https://reactome.org/ContentService/data/participants/R-HSA-191273`

```json
[
  {
    "peDbId": 196423,
    "displayName": "DHCR7 [endoplasmic reticulum membrane]",
    "schemaClass": "EntityWithAccessionedSequence",
    "refEntities": [{
      "dbId": 53616,
      "stId": "uniprot:Q9UBM7",
      "identifier": "Q9UBM7",
      "displayName": "UniProt:Q9UBM7 DHCR7"
    }]
  },
  {
    "peDbId": 193384,
    "displayName": "CHOL [endoplasmic reticulum membrane]",
    "schemaClass": "SimpleEntity",
    "refEntities": [{
      "stId": "chebi:16113",
      "identifier": "16113",
      "displayName": "cholesterol [ChEBI:16113]"
    }]
  },
  {
    "peDbId": 29364,
    "displayName": "NADPH [cytosol]",
    "schemaClass": "SimpleEntity",
    "refEntities": [{
      "stId": "chebi:57783",
      "identifier": "57783",
      "displayName": "NADPH(4-) [ChEBI:57783]"
    }]
  },
  {
    "peDbId": 194634,
    "displayName": "LAN [endoplasmic reticulum membrane]",
    "schemaClass": "SimpleEntity",
    "refEntities": [{
      "stId": "chebi:16521",
      "identifier": "16521",
      "displayName": "lanosterol [ChEBI:16521]"
    }]
  }
]
```

---

## 2. WikiPathways Sample Data

### Pathway Info: COVID-19 Pathway (WP4846)

**API Call:** `GET https://webservice.wikipathways.org/getPathwayInfo?pwId=WP4846&format=json`

```json
{
  "pathwayInfo": {
    "id": "WP4846",
    "url": "https://classic.wikipathways.org/index.php/Pathway:WP4846",
    "name": "SARS-CoV-2 and COVID-19 pathway",
    "species": "Homo sapiens",
    "revision": "140186"
  }
}
```

### Pathway Search: Cholesterol Pathways

**API Call:** `GET https://webservice.wikipathways.org/findPathwaysByText?query=cholesterol&species=Homo%20sapiens&format=json`

```json
{
  "result": [
    {
      "score": {"0": "4.837"},
      "id": "WP5333",
      "name": "Enterocyte cholesterol metabolism",
      "species": "Homo sapiens",
      "revision": "139904"
    },
    {
      "score": {"0": "4.797"},
      "id": "WP5304",
      "name": "Cholesterol metabolism",
      "species": "Homo sapiens",
      "revision": "141099"
    },
    {
      "score": {"0": "4.655"},
      "id": "WP197",
      "name": "Cholesterol biosynthesis pathway",
      "species": "Homo sapiens",
      "revision": "141096"
    },
    {
      "score": {"0": "4.606"},
      "id": "WP5329",
      "name": "Cholesterol biosynthesis pathway in hepatocytes",
      "species": "Homo sapiens",
      "revision": "140405"
    },
    {
      "score": {"0": "4.472"},
      "id": "WP4804",
      "name": "Cholesterol biosynthesis with skeletal dysplasias",
      "species": "Homo sapiens",
      "revision": "141098"
    }
  ]
}
```

### Curation Tags: Cholesterol Biosynthesis (WP197)

**API Call:** `GET https://webservice.wikipathways.org/getCurationTags?pwId=WP197&format=json`

```json
{
  "tags": [
    {
      "name": "Curation:AnalysisCollection",
      "displayName": "Approved version",
      "pathway": {
        "id": "WP197",
        "name": "Cholesterol biosynthesis pathway",
        "species": "Homo sapiens"
      },
      "revision": "141096",
      "timeModified": "20251103121712",
      "userModified": "Eweitz"
    },
    {
      "name": "Curation:FeaturedPathway",
      "displayName": "Featured version",
      "pathway": {
        "id": "WP197",
        "name": "Cholesterol biosynthesis pathway",
        "species": "Homo sapiens"
      },
      "revision": "141096"
    }
  ]
}
```

### Organism List

**API Call:** `GET https://webservice.wikipathways.org/listOrganisms?format=json`

```json
{
  "organisms": [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Danio rerio",
    "Drosophila melanogaster",
    "Caenorhabditis elegans",
    "Saccharomyces cerevisiae",
    "Arabidopsis thaliana",
    "Bos taurus",
    "Gallus gallus",
    "Escherichia coli",
    "Zea mays"
  ]
}
```

### GPML Sample: Cholesterol Biosynthesis Pathway

**API Call:** `GET https://webservice.wikipathways.org/getPathway?pwId=WP197&format=json`

```xml
<?xml version="1.0" encoding="UTF-8"?>
<Pathway xmlns="http://pathvisio.org/GPML/2013a"
         Name="Cholesterol biosynthesis pathway"
         Organism="Homo sapiens">

  <Comment Source="WikiPathways-description">
    Cholesterol is a waxy steroid metabolite found in cell membranes...
  </Comment>

  <Graphics BoardWidth="642.628" BoardHeight="648.152" />

  <DataNode TextLabel="HMGCR" GraphId="d8b" Type="GeneProduct">
    <Graphics CenterX="218.39" CenterY="190.95"
              Width="66.667" Height="20.0" ZOrder="32768"
              FontSize="12" Valign="Middle" />
    <Xref Database="Entrez Gene" ID="3156" />
  </DataNode>

  <DataNode TextLabel="Cholesterol" GraphId="e9f" Type="Metabolite">
    <Graphics CenterX="497.109" CenterY="80.142"
              Width="66.667" Height="20.0" ZOrder="32768"
              FontSize="12" Valign="Middle" Color="0000ff" />
    <Xref Database="CAS" ID="57-88-5" />
  </DataNode>

  <DataNode TextLabel="HMG-CoA" GraphId="e89f4" Type="Metabolite">
    <Graphics CenterX="300.916" CenterY="155.899"
              Width="66.667" Height="20.0" ZOrder="32768"
              FontSize="12" Valign="Middle" Color="0000ff" />
    <Xref Database="HMDB" ID="HMDB0001375" />
  </DataNode>

  <DataNode TextLabel="Mevalonic acid" GraphId="a6a21" Type="Metabolite">
    <Graphics CenterX="300.916" CenterY="227.249"
              Width="90.0" Height="20.0" ZOrder="32768"
              FontSize="12" Valign="Middle" Color="0000ff" />
    <Xref Database="CAS" ID="150-97-0" />
  </DataNode>

  <DataNode TextLabel="Squalene" GraphId="f02" Type="Metabolite">
    <Graphics CenterX="497.109" CenterY="478.127"
              Width="66.667" Height="20.0" ZOrder="32768"
              FontSize="12" Valign="Middle" Color="0000ff" />
    <Xref Database="CAS" ID="111-02-4" />
  </DataNode>

  <Interaction GraphId="a94b2">
    <Graphics ZOrder="12288" LineThickness="1.0">
      <Point X="300.916" Y="90.142" GraphRef="fb3" RelX="0.0" RelY="1.0" />
      <Point X="300.916" Y="145.899" GraphRef="e89f4" RelX="0.0" RelY="-1.0"
             ArrowHead="mim-conversion" />
    </Graphics>
    <Xref Database="" ID="" />
  </Interaction>

</Pathway>
```

---

## 3. DisGeNET Sample Data

### Gene-Disease Association

Based on DisGeNET API documentation and schema:

```json
{
  "geneId": 3156,
  "geneSymbol": "HMGCR",
  "geneName": "3-hydroxy-3-methylglutaryl-CoA reductase",
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "score": 0.82,
  "EI": 1.0,
  "DSI": 0.38,
  "DPI": 0.71,
  "associationType": "GeneticVariation",
  "source": "CTD_human",
  "NofPmids": 245,
  "pmids": [12345678, 23456789, 34567890]
}
```

### Variant-Disease Association

```json
{
  "snpId": "rs3846662",
  "chromosome": "5",
  "position": 75360714,
  "reference": "A",
  "alternative": "G",
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "score": 0.65,
  "EI": 0.92,
  "consequence": "intron_variant",
  "geneSymbol": "HMGCR",
  "source": "GWAS Catalog",
  "gnomAD_AF": 0.42,
  "NofPmids": 18
}
```

### Gene Attributes

```json
{
  "geneId": 3156,
  "geneSymbol": "HMGCR",
  "geneName": "3-hydroxy-3-methylglutaryl-CoA reductase",
  "uniprotId": ["P04035"],
  "proteinClass": "Enzyme",
  "pLI": 0.12,
  "DSI": 0.38,
  "DPI": 0.71,
  "NofDiseases": 156
}
```

### Disease Attributes

```json
{
  "diseaseId": "C0020443",
  "diseaseName": "Hypercholesterolemia",
  "diseaseType": "Disease",
  "diseaseClass": ["Metabolic Diseases", "Nutritional and Metabolic Diseases"],
  "diseaseSemanticType": "Disease or Syndrome",
  "DOID": "DOID:1168",
  "HPO": ["HP:0003124"],
  "NofGenes": 312,
  "NofVariants": 1847
}
```

---

## 4. Cross-Database Relationships

### Cholesterol Biosynthesis Across Databases

| Entity | Reactome | WikiPathways | DisGeNET |
|--------|----------|--------------|----------|
| **Pathway** | R-HSA-191273 | WP197 | - |
| **HMGCR Gene** | UniProt:P04035 | Entrez:3156 | geneId:3156 |
| **Cholesterol** | ChEBI:16113 | CAS:57-88-5 | - |
| **Hypercholesterolemia** | - | - | C0020443 |

### Identifier Cross-References

| Database | HMGCR IDs |
|----------|-----------|
| NCBI Gene | 3156 |
| UniProt | P04035 |
| Ensembl | ENSG00000113161 |
| HGNC | HGNC:5006 |
| RefSeq | NM_000859 |

| Database | Cholesterol IDs |
|----------|-----------------|
| ChEBI | CHEBI:16113 |
| HMDB | HMDB0000067 |
| KEGG | C00187 |
| PubChem | 5997 |
| CAS | 57-88-5 |

---

## 5. Database Statistics Summary

| Database | Version | Pathways | Genes | Molecules | Diseases |
|----------|---------|----------|-------|-----------|----------|
| Reactome | 95 | 2,712 | 11,196 | 1,925 | - |
| WikiPathways | 2026-01 | 3,100+ | - | - | - |
| DisGeNET | 7.0 | - | 17,549 | - | 24,166 |
