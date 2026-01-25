---
id: schema-pubmed-central
title: "PubMed Central (PMC) Schema Documentation"
type: schema
parent: README.md
last_updated: 2026-01-24
status: final
tags: [schema, database, literature, full-text]
---

# PubMed Central (PMC) - Schema Documentation

## TL;DR

PubMed Central stores full-text articles in JATS XML (Journal Article Tag Suite) format, enabling structured access to article sections, figures, tables, and supplementary materials.

## JATS XML Structure

### Document Outline

```xml
<?xml version="1.0"?>
<!DOCTYPE article PUBLIC "-//NLM//DTD JATS (Z39.96) Journal Publishing DTD v1.3">
<article article-type="research-article">
  <front>
    <!-- Metadata -->
  </front>
  <body>
    <!-- Main content -->
  </body>
  <back>
    <!-- References, acknowledgments -->
  </back>
  <floats-group>
    <!-- Figures, tables -->
  </floats-group>
</article>
```

### Front Matter (Metadata)

```xml
<front>
  <journal-meta>
    <journal-id journal-id-type="nlm-ta">Nat Commun</journal-id>
    <journal-title-group>
      <journal-title>Nature Communications</journal-title>
    </journal-title-group>
    <issn pub-type="epub">2041-1723</issn>
    <publisher>
      <publisher-name>Nature Publishing Group</publisher-name>
    </publisher>
  </journal-meta>

  <article-meta>
    <article-id pub-id-type="pmid">12345678</article-id>
    <article-id pub-id-type="pmc">PMC1234567</article-id>
    <article-id pub-id-type="doi">10.1038/ncomms12345</article-id>

    <article-categories>
      <subj-group subj-group-type="heading">
        <subject>Article</subject>
      </subj-group>
    </article-categories>

    <title-group>
      <article-title>Article Title Here</article-title>
    </title-group>

    <contrib-group>
      <contrib contrib-type="author">
        <name>
          <surname>Smith</surname>
          <given-names>John A.</given-names>
        </name>
        <xref ref-type="aff" rid="aff1">1</xref>
        <contrib-id contrib-id-type="orcid">0000-0002-1825-0097</contrib-id>
      </contrib>
    </contrib-group>

    <aff id="aff1">
      <institution>Harvard University</institution>
      <addr-line>Cambridge, MA, USA</addr-line>
    </aff>

    <pub-date pub-type="epub">
      <day>15</day>
      <month>01</month>
      <year>2024</year>
    </pub-date>

    <abstract>
      <p>Abstract text here...</p>
    </abstract>

    <kwd-group>
      <kwd>keyword1</kwd>
      <kwd>keyword2</kwd>
    </kwd-group>
  </article-meta>
</front>
```

### Body Content

```xml
<body>
  <sec id="sec1">
    <title>Introduction</title>
    <p>Introduction text with <xref ref-type="bibr" rid="R1">references</xref>.</p>
    <p>Text mentioning <xref ref-type="fig" rid="fig1">Figure 1</xref>.</p>
  </sec>

  <sec id="sec2">
    <title>Results</title>
    <sec id="sec2.1">
      <title>Subsection</title>
      <p>Content with <italic>italic</italic> and <bold>bold</bold> text.</p>
      <p>Gene names: <named-content content-type="gene">TP53</named-content></p>
    </sec>
  </sec>

  <sec id="sec3">
    <title>Methods</title>
    <p>Methods description...</p>
  </sec>
</body>
```

### Figures and Tables

```xml
<fig id="fig1">
  <label>Figure 1</label>
  <caption>
    <title>Figure title</title>
    <p>Detailed caption text.</p>
  </caption>
  <graphic xlink:href="ncomms12345-f1.jpg"/>
</fig>

<table-wrap id="tbl1">
  <label>Table 1</label>
  <caption>
    <title>Table title</title>
  </caption>
  <table>
    <thead>
      <tr>
        <th>Column 1</th>
        <th>Column 2</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>Data 1</td>
        <td>Data 2</td>
      </tr>
    </tbody>
  </table>
</table-wrap>
```

### Back Matter (References)

```xml
<back>
  <ack>
    <title>Acknowledgements</title>
    <p>We thank...</p>
  </ack>

  <fn-group>
    <fn fn-type="financial-disclosure">
      <p>Funding: NIH grant R01-GM123456</p>
    </fn>
  </fn-group>

  <ref-list>
    <title>References</title>
    <ref id="R1">
      <label>1</label>
      <element-citation publication-type="journal">
        <person-group person-group-type="author">
          <name>
            <surname>Doe</surname>
            <given-names>J</given-names>
          </name>
        </person-group>
        <article-title>Reference title</article-title>
        <source>Journal Name</source>
        <year>2023</year>
        <volume>10</volume>
        <fpage>100</fpage>
        <lpage>110</lpage>
        <pub-id pub-id-type="pmid">23456789</pub-id>
        <pub-id pub-id-type="doi">10.1000/example</pub-id>
      </element-citation>
    </ref>
  </ref-list>

  <sec sec-type="supplementary-material">
    <title>Supplementary Information</title>
    <supplementary-material id="supp1">
      <label>Supplementary Table 1</label>
      <media xlink:href="supp_table1.xlsx"/>
    </supplementary-material>
  </sec>
</back>
```

## Article Identifiers

| ID Type | Attribute | Example |
|---------|-----------|---------|
| PMID | `pub-id-type="pmid"` | 12345678 |
| PMCID | `pub-id-type="pmc"` | PMC1234567 |
| DOI | `pub-id-type="doi"` | 10.1038/example |
| Publisher ID | `pub-id-type="publisher-id"` | ncomms12345 |
| Manuscript ID | `pub-id-type="manuscript"` | NIHMS123456 |

## Article Types

| Type | Description |
|------|-------------|
| research-article | Original research |
| review-article | Review paper |
| case-report | Clinical case report |
| letter | Letter to editor |
| editorial | Editorial content |
| correction | Correction/erratum |
| retraction | Retraction notice |

## Cross-Reference Types

| ref-type | Target |
|----------|--------|
| bibr | Bibliography reference |
| fig | Figure |
| table | Table |
| supplementary-material | Supplement |
| aff | Affiliation |
| sec | Section |

## Named Content Types

| content-type | Description |
|--------------|-------------|
| gene | Gene names |
| species | Organism names |
| compound | Chemical compounds |
| disease | Disease terms |

## OA Subset License Tags

```xml
<license license-type="open-access">
  <ali:license_ref>https://creativecommons.org/licenses/by/4.0/</ali:license_ref>
  <license-p>This is an open access article...</license-p>
</license>
```

## File List Format (OA Subset)

```
oa_comm_use_file_list.txt:
PMC1234567<TAB>J Biol Chem<TAB>file.tar.gz<TAB>CC BY
PMC1234568<TAB>Nature<TAB>file.tar.gz<TAB>CC BY-NC
```

## See Also

- [Download Documentation](./download.md)
- [JATS Standard](https://jats.nlm.nih.gov/)
