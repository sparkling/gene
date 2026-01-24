---
id: schema-pmc-id-converter
title: "PMC ID Converter Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, database, identifier-mapping, literature]
---

# PMC ID Converter - Schema Documentation

## TL;DR

The PMC ID Converter service converts between publication identifier types (PMID, PMCID, DOI) and returns additional metadata about articles. Responses are available in JSON, XML, and CSV formats.

## Response Formats

### JSON Response

```json
{
  "status": "ok",
  "responseDate": "2024-01-15 10:30:00",
  "request": {
    "ids": "12345678",
    "idtype": "pmid",
    "format": "json",
    "versions": "yes"
  },
  "records": [
    {
      "pmid": "12345678",
      "pmcid": "PMC1234567",
      "doi": "10.1038/example",
      "live": true,
      "status": "current",
      "errata": null,
      "versions": [
        {
          "pmcid": "PMC1234567.1",
          "current": "true"
        }
      ]
    }
  ]
}
```

### XML Response

```xml
<?xml version="1.0"?>
<!DOCTYPE pmcids PUBLIC "-//NLM//DTD PMC ID Converter Output//EN" "pmcids.dtd">
<pmcids status="ok">
  <request>
    <ids>12345678</ids>
    <idtype>pmid</idtype>
    <format>xml</format>
    <versions>yes</versions>
  </request>
  <record requested-id="12345678" pmcid="PMC1234567" pmid="12345678" doi="10.1038/example">
    <versions>
      <version pmcid="PMC1234567.1" current="true" />
    </versions>
  </record>
</pmcids>
```

### CSV Response

```csv
PMID,PMCID,DOI,Status
12345678,PMC1234567,10.1038/example,current
23456789,PMC2345678,10.1000/example2,current
```

## Record Fields

| Field | Type | Description |
|-------|------|-------------|
| pmid | string | PubMed identifier |
| pmcid | string | PMC identifier (with PMC prefix) |
| doi | string | Digital Object Identifier |
| mid | string | Manuscript ID (NIHMS format) |
| live | boolean | Article is publicly accessible |
| status | string | Article status |
| errata | object | Correction information (if any) |
| versions | array | Version history |
| release-date | string | Public release date |

## Status Values

| Status | Description |
|--------|-------------|
| current | Live, current version |
| retracted | Article has been retracted |
| removed | Article removed from PMC |
| preprint | Preprint version |

## Identifier Patterns

| Type | Pattern | Example |
|------|---------|---------|
| PMID | `[0-9]+` | 12345678 |
| PMCID | `PMC[0-9]+` | PMC1234567 |
| DOI | `10\.[0-9]+/.+` | 10.1038/example |
| MID | `NIHMS[0-9]+` | NIHMS123456 |

## Version Object

```json
{
  "pmcid": "PMC1234567.1",
  "current": "true",
  "release-date": "2024-01-15"
}
```

Versions are indicated by a decimal suffix (e.g., PMC1234567.1, PMC1234567.2).

## Errata Object

When an article has corrections:

```json
{
  "errata": {
    "type": "correction",
    "pmid": "34567890",
    "pmcid": "PMC3456789",
    "doi": "10.1038/correction"
  }
}
```

## Batch Response

```json
{
  "status": "ok",
  "records": [
    {
      "pmid": "12345678",
      "pmcid": "PMC1234567",
      "doi": "10.1038/example1"
    },
    {
      "pmid": "23456789",
      "pmcid": "PMC2345678",
      "doi": "10.1038/example2"
    },
    {
      "pmid": "34567890",
      "errmsg": "No PMC record found"
    }
  ]
}
```

## Error Responses

### Invalid ID

```json
{
  "status": "ok",
  "records": [
    {
      "requested-id": "invalid123",
      "errmsg": "No PMID found"
    }
  ]
}
```

### No PMC Record

```json
{
  "records": [
    {
      "pmid": "12345678",
      "errmsg": "No PMC record found"
    }
  ]
}
```

### Rate Limit Exceeded

```json
{
  "status": "error",
  "message": "Too many requests. Please try again later."
}
```

## Request Parameters

| Parameter | Values | Description |
|-----------|--------|-------------|
| ids | comma-separated | Input identifiers (max 200) |
| idtype | pmid, pmcid, doi, mid | Input ID type (auto-detect if omitted) |
| format | json, xml, csv | Response format |
| versions | yes, no | Include version history |
| tool | string | Your application name |
| email | string | Contact email (recommended) |

## Release Date Format

Dates are returned in ISO format: `YYYY-MM-DD`

## See Also

- [Download Documentation](./download.md)
- [PMC ID Converter API](https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/)
