---
id: schema-pubmed
title: "PubMed Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, literature, biomedical]
---

# PubMed - Schema Documentation

## Entity Relationship Overview

```
┌───────────────┐              ┌───────────────┐
│    Journal    │<─────────────│    Article    │
│   (ISSN)      │              │    (PMID)     │
└───────────────┘              └───────┬───────┘
                                       │
       ┌───────────────┬───────────────┼───────────────┬───────────────┐
       │               │               │               │               │
       v               v               v               v               v
┌───────────┐   ┌───────────┐   ┌───────────┐   ┌───────────┐   ┌───────────┐
│  Author   │   │  Abstract │   │   MeSH    │   │  Keyword  │   │ Reference │
│ (ORCID)   │   │  (Text)   │   │  (Terms)  │   │  (Author) │   │  (PMID)   │
└───────────┘   └───────────┘   └─────┬─────┘   └───────────┘   └───────────┘
      │                               │
      v                               v
┌───────────┐                 ┌───────────────┐
│Affiliation│                 │  Qualifier    │
│  (Org)    │                 │ (Subheading)  │
└───────────┘                 └───────────────┘

Relationships:
- Journal (1) ----< (N) Article: One journal contains many articles
- Article (1) ----< (N) Author: One article has many authors
- Article (1) ----< (N) MeSH Heading: One article has many MeSH terms
- Article (1) ----< (N) Reference: One article cites many references
- MeSH Heading (1) ----< (N) Qualifier: One heading has many qualifiers
- Author (N) >----< (M) Affiliation: Authors linked to institutions
```

---

## TL;DR

PubMed provides bibliographic records in multiple XML formats including PubMed XML, MEDLINE format, and JSON. Each record contains citation metadata, abstract text, MeSH terms, and publication details.

## PubMed Article Set (XML)

### Document Structure

```xml
<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation Status="MEDLINE" Owner="NLM">
      <!-- Article metadata -->
    </MedlineCitation>
    <PubmedData>
      <!-- Publication history, IDs -->
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
```

### MedlineCitation Element

```xml
<MedlineCitation Status="MEDLINE" Owner="NLM">
  <PMID Version="1">12345678</PMID>
  <DateCompleted>
    <Year>2024</Year>
    <Month>01</Month>
    <Day>15</Day>
  </DateCompleted>
  <DateRevised>
    <Year>2024</Year>
    <Month>01</Month>
    <Day>20</Day>
  </DateRevised>

  <Article PubModel="Print-Electronic">
    <Journal>
      <ISSN IssnType="Electronic">1476-4687</ISSN>
      <JournalIssue CitedMedium="Internet">
        <Volume>615</Volume>
        <Issue>7950</Issue>
        <PubDate>
          <Year>2024</Year>
          <Month>Jan</Month>
        </PubDate>
      </JournalIssue>
      <Title>Nature</Title>
      <ISOAbbreviation>Nature</ISOAbbreviation>
    </Journal>

    <ArticleTitle>Article Title Here</ArticleTitle>

    <Pagination>
      <StartPage>100</StartPage>
      <EndPage>110</EndPage>
      <MedlinePgn>100-110</MedlinePgn>
    </Pagination>

    <ELocationID EIdType="doi" ValidYN="Y">10.1038/example</ELocationID>
    <ELocationID EIdType="pii" ValidYN="Y">s41586-023-12345-6</ELocationID>

    <Abstract>
      <AbstractText>Full abstract text here...</AbstractText>
    </Abstract>

    <AuthorList CompleteYN="Y">
      <Author ValidYN="Y">
        <LastName>Smith</LastName>
        <ForeName>John A</ForeName>
        <Initials>JA</Initials>
        <Identifier Source="ORCID">0000-0002-1825-0097</Identifier>
        <AffiliationInfo>
          <Affiliation>Department of Biology, Harvard University, Cambridge, MA, USA.</Affiliation>
        </AffiliationInfo>
      </Author>
    </AuthorList>

    <Language>eng</Language>

    <PublicationTypeList>
      <PublicationType UI="D016428">Journal Article</PublicationType>
      <PublicationType UI="D013485">Research Support, Non-U.S. Gov't</PublicationType>
    </PublicationTypeList>

    <ArticleDate DateType="Electronic">
      <Year>2024</Year>
      <Month>01</Month>
      <Day>10</Day>
    </ArticleDate>
  </Article>

  <MeshHeadingList>
    <MeshHeading>
      <DescriptorName UI="D006801" MajorTopicYN="N">Humans</DescriptorName>
    </MeshHeading>
    <MeshHeading>
      <DescriptorName UI="D016159" MajorTopicYN="Y">Tumor Suppressor Protein p53</DescriptorName>
      <QualifierName UI="Q000235" MajorTopicYN="N">genetics</QualifierName>
      <QualifierName UI="Q000378" MajorTopicYN="N">metabolism</QualifierName>
    </MeshHeading>
  </MeshHeadingList>

  <KeywordList Owner="NOTNLM">
    <Keyword MajorTopicYN="N">TP53</Keyword>
    <Keyword MajorTopicYN="N">cancer</Keyword>
  </KeywordList>
</MedlineCitation>
```

### PubmedData Element

```xml
<PubmedData>
  <History>
    <PubMedPubDate PubStatus="received">
      <Year>2023</Year>
      <Month>06</Month>
      <Day>15</Day>
    </PubMedPubDate>
    <PubMedPubDate PubStatus="accepted">
      <Year>2023</Year>
      <Month>12</Month>
      <Day>01</Day>
    </PubMedPubDate>
    <PubMedPubDate PubStatus="pubmed">
      <Year>2024</Year>
      <Month>01</Month>
      <Day>11</Day>
    </PubMedPubDate>
    <PubMedPubDate PubStatus="medline">
      <Year>2024</Year>
      <Month>01</Month>
      <Day>16</Day>
    </PubMedPubDate>
    <PubMedPubDate PubStatus="entrez">
      <Year>2024</Year>
      <Month>01</Month>
      <Day>10</Day>
    </PubMedPubDate>
  </History>

  <PublicationStatus>ppublish</PublicationStatus>

  <ArticleIdList>
    <ArticleId IdType="pubmed">12345678</ArticleId>
    <ArticleId IdType="pmc">PMC1234567</ArticleId>
    <ArticleId IdType="doi">10.1038/example</ArticleId>
    <ArticleId IdType="pii">s41586-023-12345-6</ArticleId>
  </ArticleIdList>

  <ReferenceList>
    <Reference>
      <Citation>Nature. 2020;580(7803):360-366</Citation>
      <ArticleIdList>
        <ArticleId IdType="pubmed">32322060</ArticleId>
      </ArticleIdList>
    </Reference>
  </ReferenceList>
</PubmedData>
```

## Citation Status Values

| Status | Description |
|--------|-------------|
| MEDLINE | Fully indexed with MeSH |
| PubMed-not-MEDLINE | Not indexed (e-ahead, etc.) |
| In-Process | Being processed |
| Publisher | Direct from publisher |
| Preprint | Preprint record |

## Publication Types

| UI | Type |
|----|------|
| D016428 | Journal Article |
| D016454 | Review |
| D017418 | Meta-Analysis |
| D016449 | Randomized Controlled Trial |
| D016422 | Letter |
| D016420 | Comment |
| D016421 | Editorial |
| D002363 | Case Reports |

## MeSH Term Structure

| Attribute | Description |
|-----------|-------------|
| UI | Unique identifier (D######) |
| MajorTopicYN | Y = major topic |
| QualifierName | Subheading (optional) |

## JSON Format (E-utilities)

```json
{
  "header": {
    "type": "esummary",
    "version": "0.3"
  },
  "result": {
    "uids": ["12345678"],
    "12345678": {
      "uid": "12345678",
      "pubdate": "2024 Jan",
      "epubdate": "2024 Jan 10",
      "source": "Nature",
      "authors": [
        {
          "name": "Smith JA",
          "authtype": "Author",
          "clusterid": ""
        }
      ],
      "lastauthor": "Smith JA",
      "title": "Article Title Here",
      "sortitle": "article title here",
      "volume": "615",
      "issue": "7950",
      "pages": "100-110",
      "lang": ["eng"],
      "nlmuniqueid": "0410462",
      "issn": "0028-0836",
      "essn": "1476-4687",
      "pubtype": ["Journal Article"],
      "recordstatus": "PubMed - indexed for MEDLINE",
      "pubstatus": "256",
      "articleids": [
        {"idtype": "pubmed", "value": "12345678"},
        {"idtype": "pmc", "value": "PMC1234567"},
        {"idtype": "doi", "value": "10.1038/example"}
      ],
      "fulljournalname": "Nature",
      "sortpubdate": "2024/01/01 00:00"
    }
  }
}
```

## Structured Abstract

```xml
<Abstract>
  <AbstractText Label="BACKGROUND" NlmCategory="BACKGROUND">
    Background text...
  </AbstractText>
  <AbstractText Label="METHODS" NlmCategory="METHODS">
    Methods text...
  </AbstractText>
  <AbstractText Label="RESULTS" NlmCategory="RESULTS">
    Results text...
  </AbstractText>
  <AbstractText Label="CONCLUSIONS" NlmCategory="CONCLUSIONS">
    Conclusions text...
  </AbstractText>
</Abstract>
```

## Grant Information

```xml
<GrantList CompleteYN="Y">
  <Grant>
    <GrantID>R01 GM123456</GrantID>
    <Acronym>GM</Acronym>
    <Agency>NIGMS NIH HHS</Agency>
    <Country>United States</Country>
  </Grant>
</GrantList>
```

## See Also

- [Download Documentation](./download.md)
- [PubMed DTD](https://dtd.nlm.nih.gov/)
