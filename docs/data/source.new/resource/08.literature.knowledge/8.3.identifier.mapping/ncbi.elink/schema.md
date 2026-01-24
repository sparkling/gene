---
id: schema-ncbi-elink
title: "NCBI E-Link Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, database, identifier-mapping, ncbi]
---

# NCBI E-Link - Schema Documentation

## TL;DR

NCBI E-Link (ELink) returns XML responses describing links between records across NCBI databases. It supports multiple link types and can return related records, external links, and database-specific connections.

## Response Structure

### Basic Link Response

```xml
<?xml version="1.0"?>
<!DOCTYPE eLinkResult PUBLIC "-//NLM//DTD elink 20101123//EN" "https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd">
<eLinkResult>
  <LinkSet>
    <DbFrom>pubmed</DbFrom>
    <IdList>
      <Id>12345678</Id>
    </IdList>
    <LinkSetDb>
      <DbTo>gene</DbTo>
      <LinkName>pubmed_gene</LinkName>
      <Link>
        <Id>7157</Id>
        <Score>1000</Score>
      </Link>
      <Link>
        <Id>672</Id>
        <Score>800</Score>
      </Link>
    </LinkSetDb>
  </LinkSet>
</eLinkResult>
```

### JSON Format

```json
{
  "header": {
    "type": "elink",
    "version": "0.3"
  },
  "linksets": [
    {
      "dbfrom": "pubmed",
      "ids": ["12345678"],
      "linksetdbs": [
        {
          "dbto": "gene",
          "linkname": "pubmed_gene",
          "links": ["7157", "672"]
        }
      ]
    }
  ]
}
```

## Link Types

### cmd=neighbor (Related Records)

```xml
<LinkSetDb>
  <DbTo>pubmed</DbTo>
  <LinkName>pubmed_pubmed</LinkName>
  <Link>
    <Id>23456789</Id>
    <Score>50234</Score>
  </Link>
</LinkSetDb>
```

### cmd=neighbor_score (With Similarity Scores)

```xml
<LinkSetDb>
  <DbTo>pubmed</DbTo>
  <LinkName>pubmed_pubmed</LinkName>
  <Link>
    <Id>23456789</Id>
    <Score>95</Score>
  </Link>
</LinkSetDb>
```

### cmd=acheck (Available Links)

```xml
<LinkSet>
  <DbFrom>pubmed</DbFrom>
  <IdCheckList>
    <Id>
      <HasLinkOut>Y</HasLinkOut>
      <HasNeighbor>Y</HasNeighbor>
      <value>12345678</value>
    </Id>
  </IdCheckList>
</LinkSet>
```

### cmd=llinks (External Links)

```xml
<LinkSet>
  <IdUrlList>
    <IdUrlSet>
      <Id>12345678</Id>
      <ObjUrl>
        <Url>https://www.example.com/article/12345678</Url>
        <LinkName>Full Text</LinkName>
        <SubjectType>full-text</SubjectType>
        <Category>Full Text Sources</Category>
        <Attribute>free resource</Attribute>
        <Provider>
          <Name>Example Publisher</Name>
          <NameAbbr>EP</NameAbbr>
          <Id>1234</Id>
          <Url>https://www.example.com</Url>
        </Provider>
      </ObjUrl>
    </IdUrlSet>
  </IdUrlList>
</LinkSet>
```

## Common Link Names

### PubMed Links

| Link Name | From | To | Description |
|-----------|------|-----|-------------|
| pubmed_pubmed | pubmed | pubmed | Related articles |
| pubmed_gene | pubmed | gene | Associated genes |
| pubmed_protein | pubmed | protein | Associated proteins |
| pubmed_pmc | pubmed | pmc | Full text in PMC |
| pubmed_mesh | pubmed | mesh | MeSH terms |
| pubmed_structure | pubmed | structure | 3D structures |

### Gene Links

| Link Name | From | To | Description |
|-----------|------|-----|-------------|
| gene_pubmed | gene | pubmed | Gene literature |
| gene_protein | gene | protein | Gene products |
| gene_nuccore | gene | nuccore | Nucleotide sequences |
| gene_snp | gene | snp | SNP variants |
| gene_homologene | gene | homologene | Homologs |
| gene_clinvar | gene | clinvar | Clinical variants |

### Protein Links

| Link Name | From | To | Description |
|-----------|------|-----|-------------|
| protein_pubmed | protein | pubmed | Protein literature |
| protein_gene | protein | gene | Encoding gene |
| protein_structure | protein | structure | 3D structures |
| protein_cdd | protein | cdd | Conserved domains |
| protein_biosystems | protein | biosystems | Pathway involvement |

## Database Codes

| Code | Database | Description |
|------|----------|-------------|
| pubmed | PubMed | Literature |
| pmc | PubMed Central | Full text |
| gene | Gene | Gene records |
| protein | Protein | Protein sequences |
| nuccore | Nucleotide | Nucleotide sequences |
| structure | Structure | 3D structures |
| snp | dbSNP | Variants |
| clinvar | ClinVar | Clinical variants |
| mesh | MeSH | Vocabulary |
| taxonomy | Taxonomy | Species |
| biosystems | BioSystems | Pathways |
| omim | OMIM | Genetic disorders |
| gds | GEO DataSets | Expression |
| homologene | HomoloGene | Homologs |

## Batch Response

```xml
<eLinkResult>
  <LinkSet>
    <DbFrom>pubmed</DbFrom>
    <IdList>
      <Id>12345678</Id>
    </IdList>
    <LinkSetDb>
      <DbTo>gene</DbTo>
      <LinkName>pubmed_gene</LinkName>
      <Link><Id>7157</Id></Link>
    </LinkSetDb>
  </LinkSet>
  <LinkSet>
    <DbFrom>pubmed</DbFrom>
    <IdList>
      <Id>23456789</Id>
    </IdList>
    <LinkSetDb>
      <DbTo>gene</DbTo>
      <LinkName>pubmed_gene</LinkName>
      <Link><Id>672</Id></Link>
    </LinkSetDb>
  </LinkSet>
</eLinkResult>
```

## History Server Response

When using `cmd=neighbor_history`:

```xml
<LinkSet>
  <DbFrom>pubmed</DbFrom>
  <LinkSetDbHistory>
    <DbTo>pubmed</DbTo>
    <LinkName>pubmed_pubmed</LinkName>
    <QueryKey>1</QueryKey>
  </LinkSetDbHistory>
</LinkSet>
```

## Error Response

```xml
<eLinkResult>
  <ERROR>Invalid db name</ERROR>
</eLinkResult>
```

Or for individual IDs:

```xml
<LinkSet>
  <DbFrom>pubmed</DbFrom>
  <IdList>
    <Id>99999999999</Id>
  </IdList>
  <ERROR>Id not found</ERROR>
</LinkSet>
```

## See Also

- [Download Documentation](./download.md)
- [E-utilities Documentation](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
