# Research Paper Data Structures for Genetics Knowledge Base

**Last Updated:** January 2026
**Purpose:** Data formats, parsing strategies, and embedding approaches for biomedical research papers
**Target System:** RuVector (vector database with Cypher graph queries)

---

## Overview

Research papers are critical for the genetics knowledge base, providing:
1. **Evidence** - Citations for SNP-phenotype associations
2. **Context** - Detailed mechanisms and study findings
3. **Discovery** - New relationships via semantic search and RAG
4. **Authority** - Academic credibility for recommendations

This document covers data structures from acquisition through storage in RuVector.

---

## 1. PubMed XML Format

### 1.1 DTD Structure

PubMed uses the `pubmed_230101.dtd` Document Type Definition (updated annually). The core structure:

```xml
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMed 230101//EN"
    "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_230101.dtd">
```

### 1.2 Key XML Elements

```xml
<?xml version="1.0" encoding="UTF-8"?>
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation Status="MEDLINE" Owner="NLM">
      <!-- Primary identifier -->
      <PMID Version="1">12345678</PMID>

      <DateCompleted>
        <Year>2024</Year>
        <Month>03</Month>
        <Day>15</Day>
      </DateCompleted>

      <Article PubModel="Print-Electronic">
        <!-- Journal information -->
        <Journal>
          <ISSN IssnType="Electronic">1476-4687</ISSN>
          <JournalIssue CitedMedium="Internet">
            <Volume>620</Volume>
            <Issue>7975</Issue>
            <PubDate>
              <Year>2024</Year>
              <Month>Aug</Month>
            </PubDate>
          </JournalIssue>
          <Title>Nature</Title>
          <ISOAbbreviation>Nature</ISOAbbreviation>
        </Journal>

        <!-- Title -->
        <ArticleTitle>MTHFR C677T polymorphism and folate metabolism:
          A comprehensive genetic analysis</ArticleTitle>

        <!-- Pagination -->
        <Pagination>
          <MedlinePgn>234-241</MedlinePgn>
        </Pagination>

        <!-- Abstract with structured labels -->
        <Abstract>
          <AbstractText Label="BACKGROUND">
            The MTHFR C677T polymorphism affects folate metabolism...
          </AbstractText>
          <AbstractText Label="METHODS">
            We conducted a genome-wide association study...
          </AbstractText>
          <AbstractText Label="RESULTS">
            Carriers of the T allele showed 30-70% reduced enzyme activity...
          </AbstractText>
          <AbstractText Label="CONCLUSIONS">
            This variant has significant implications for...
          </AbstractText>
        </Abstract>

        <!-- Author list -->
        <AuthorList CompleteYN="Y">
          <Author ValidYN="Y">
            <LastName>Smith</LastName>
            <ForeName>John A</ForeName>
            <Initials>JA</Initials>
            <Identifier Source="ORCID">0000-0001-2345-6789</Identifier>
            <AffiliationInfo>
              <Affiliation>Department of Genetics, Harvard Medical School,
                Boston, MA, USA</Affiliation>
            </AffiliationInfo>
          </Author>
          <Author ValidYN="Y" EqualContrib="Y">
            <LastName>Chen</LastName>
            <ForeName>Wei</ForeName>
            <Initials>W</Initials>
            <AffiliationInfo>
              <Affiliation>Beijing Institute of Genomics, Chinese Academy
                of Sciences, Beijing, China</Affiliation>
            </AffiliationInfo>
          </Author>
        </AuthorList>

        <!-- Language -->
        <Language>eng</Language>

        <!-- Publication types -->
        <PublicationTypeList>
          <PublicationType UI="D016428">Journal Article</PublicationType>
          <PublicationType UI="D016454">Review</PublicationType>
          <PublicationType UI="D017418">Meta-Analysis</PublicationType>
        </PublicationTypeList>

        <!-- Article IDs -->
        <ArticleIdList>
          <ArticleId IdType="pubmed">12345678</ArticleId>
          <ArticleId IdType="doi">10.1038/s41586-024-07123-4</ArticleId>
          <ArticleId IdType="pmc">PMC12345678</ArticleId>
          <ArticleId IdType="pii">s41586-024-07123-4</ArticleId>
        </ArticleIdList>
      </Article>

      <!-- MeSH terms (controlled vocabulary) -->
      <MeshHeadingList>
        <MeshHeading>
          <DescriptorName UI="D008780" MajorTopicYN="N">
            Methylenetetrahydrofolate Reductase (NADPH2)
          </DescriptorName>
          <QualifierName UI="Q000235" MajorTopicYN="Y">genetics</QualifierName>
        </MeshHeading>
        <MeshHeading>
          <DescriptorName UI="D011110" MajorTopicYN="N">Polymorphism, Genetic</DescriptorName>
        </MeshHeading>
        <MeshHeading>
          <DescriptorName UI="D005492" MajorTopicYN="N">Folic Acid</DescriptorName>
          <QualifierName UI="Q000378" MajorTopicYN="N">metabolism</QualifierName>
        </MeshHeading>
        <MeshHeading>
          <DescriptorName UI="D006801" MajorTopicYN="N">Humans</DescriptorName>
        </MeshHeading>
      </MeshHeadingList>

      <!-- Keywords (author-supplied) -->
      <KeywordList Owner="NOTNLM">
        <Keyword MajorTopicYN="N">MTHFR</Keyword>
        <Keyword MajorTopicYN="N">folate</Keyword>
        <Keyword MajorTopicYN="N">methylation</Keyword>
        <Keyword MajorTopicYN="N">SNP</Keyword>
        <Keyword MajorTopicYN="N">rs1801133</Keyword>
      </KeywordList>

      <!-- Chemical list -->
      <ChemicalList>
        <Chemical>
          <RegistryNumber>EC 1.5.1.20</RegistryNumber>
          <NameOfSubstance UI="D008780">
            Methylenetetrahydrofolate Reductase (NADPH2)
          </NameOfSubstance>
        </Chemical>
      </ChemicalList>

    </MedlineCitation>

    <!-- PubMed-specific data -->
    <PubmedData>
      <History>
        <PubMedPubDate PubStatus="received">
          <Year>2023</Year><Month>09</Month><Day>15</Day>
        </PubMedPubDate>
        <PubMedPubDate PubStatus="accepted">
          <Year>2024</Year><Month>01</Month><Day>20</Day>
        </PubMedPubDate>
        <PubMedPubDate PubStatus="pubmed">
          <Year>2024</Year><Month>02</Month><Day>01</Day>
        </PubMedPubDate>
      </History>
      <PublicationStatus>ppublish</PublicationStatus>
      <ArticleIdList>
        <ArticleId IdType="pubmed">12345678</ArticleId>
        <ArticleId IdType="doi">10.1038/s41586-024-07123-4</ArticleId>
      </ArticleIdList>

      <!-- References (if available) -->
      <ReferenceList>
        <Reference>
          <Citation>Nature. 2020;580:123-130</Citation>
          <ArticleIdList>
            <ArticleId IdType="pubmed">32123456</ArticleId>
          </ArticleIdList>
        </Reference>
      </ReferenceList>
    </PubmedData>
  </PubmedArticle>
</PubmedArticleSet>
```

### 1.3 Parsing Considerations

**Structured vs Unstructured Abstracts:**
```typescript
// Handle both abstract formats
function parseAbstract(abstractElement: Element): string {
  const sections = abstractElement.querySelectorAll('AbstractText');

  if (sections.length > 1) {
    // Structured abstract with labels
    return Array.from(sections).map(section => {
      const label = section.getAttribute('Label') || '';
      const text = section.textContent || '';
      return label ? `${label}: ${text}` : text;
    }).join('\n\n');
  } else {
    // Unstructured abstract
    return abstractElement.textContent || '';
  }
}
```

**Author Name Variations:**
```typescript
interface PubMedAuthor {
  lastName: string;
  foreName?: string;
  initials?: string;
  collectiveName?: string;  // For group authors
  orcid?: string;
  affiliations: string[];
}

// Some articles have collective authors
// <CollectiveName>COVID-19 Host Genetics Initiative</CollectiveName>
```

**Date Handling:**
```typescript
// PubMed has multiple date types
interface PubMedDates {
  pubDate?: Date;           // Publication date
  articleDate?: Date;       // Electronic publication
  medlineDate?: string;     // Sometimes just "2024 Spring"
  received?: Date;
  accepted?: Date;
  pubmed?: Date;           // Added to PubMed date
}

// Fallback chain: PubDate > ArticleDate > MedlineDate
```

**Character Encoding:**
```typescript
// Handle special characters in titles/abstracts
const ENTITY_MAP: Record<string, string> = {
  '&lt;': '<',
  '&gt;': '>',
  '&amp;': '&',
  '&apos;': "'",
  '&quot;': '"',
  // Greek letters common in genetics
  '&alpha;': 'alpha',
  '&beta;': 'beta',
  '&gamma;': 'gamma',
};
```

---

## 2. PMC Full-Text Format (JATS XML)

### 2.1 JATS XML Schema

PMC uses the Journal Article Tag Suite (JATS) standard:

```xml
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE article PUBLIC "-//NLM//DTD JATS (Z39.96) Journal Archiving
  and Interchange DTD v1.3 20210610//EN"
  "https://jats.nlm.nih.gov/archiving/1.3/JATS-archivearticle1-3.dtd">

<article xmlns:xlink="http://www.w3.org/1999/xlink"
         article-type="research-article">

  <front>
    <!-- Article metadata (similar to PubMed) -->
    <journal-meta>
      <journal-id journal-id-type="nlm-ta">Nat Genet</journal-id>
      <journal-id journal-id-type="iso-abbrev">Nat. Genet.</journal-id>
      <journal-title-group>
        <journal-title>Nature Genetics</journal-title>
      </journal-title-group>
      <issn pub-type="ppub">1061-4036</issn>
      <issn pub-type="epub">1546-1718</issn>
      <publisher>
        <publisher-name>Nature Publishing Group</publisher-name>
      </publisher>
    </journal-meta>

    <article-meta>
      <article-id pub-id-type="pmid">12345678</article-id>
      <article-id pub-id-type="pmc">PMC12345678</article-id>
      <article-id pub-id-type="doi">10.1038/ng.1234</article-id>

      <title-group>
        <article-title>Genome-wide association study identifies
          novel MTHFR variants</article-title>
      </title-group>

      <contrib-group>
        <contrib contrib-type="author">
          <name>
            <surname>Smith</surname>
            <given-names>John A.</given-names>
          </name>
          <xref ref-type="aff" rid="aff1">1</xref>
          <xref ref-type="corresp" rid="cor1">*</xref>
        </contrib>
      </contrib-group>

      <aff id="aff1">
        <label>1</label>
        <institution>Harvard Medical School</institution>,
        <addr-line>Boston, MA</addr-line>,
        <country>USA</country>
      </aff>

      <author-notes>
        <corresp id="cor1">
          <label>*</label>Corresponding author:
          <email>jsmith@hms.harvard.edu</email>
        </corresp>
      </author-notes>

      <abstract>
        <sec>
          <title>Background</title>
          <p>The MTHFR gene encodes...</p>
        </sec>
        <sec>
          <title>Results</title>
          <p>We identified 15 novel variants...</p>
        </sec>
      </abstract>

      <kwd-group>
        <kwd>GWAS</kwd>
        <kwd>MTHFR</kwd>
        <kwd>folate metabolism</kwd>
      </kwd-group>

    </article-meta>
  </front>

  <body>
    <!-- Full article content -->
    <sec id="sec1">
      <title>Introduction</title>
      <p>Methylenetetrahydrofolate reductase (MTHFR) plays a critical
        role in folate metabolism <xref ref-type="bibr" rid="ref1">[1]</xref>.
        The C677T polymorphism (<xref ref-type="fig" rid="fig1">Fig. 1</xref>)
        results in an alanine to valine substitution...</p>

      <p>Previous genome-wide association studies (GWAS) have identified
        <xref ref-type="table" rid="tbl1">Table 1</xref> several variants
        associated with homocysteine levels...</p>
    </sec>

    <sec id="sec2">
      <title>Methods</title>

      <sec id="sec2-1">
        <title>Study Population</title>
        <p>We recruited 50,000 participants from the UK Biobank...</p>
      </sec>

      <sec id="sec2-2">
        <title>Genotyping</title>
        <p>DNA was extracted using standard protocols and genotyped
          on the Illumina Global Screening Array...</p>
      </sec>

      <sec id="sec2-3">
        <title>Statistical Analysis</title>
        <p>Association analyses were performed using PLINK v2.0
          <xref ref-type="bibr" rid="ref15">[15]</xref> with adjustment
          for population stratification...</p>

        <!-- Inline formulas -->
        <disp-formula id="eq1">
          <tex-math id="M1">
            \text{OR} = \frac{p_1 / (1-p_1)}{p_0 / (1-p_0)}
          </tex-math>
        </disp-formula>
      </sec>
    </sec>

    <sec id="sec3">
      <title>Results</title>
      <p>We identified 15 genome-wide significant variants
        (P &lt; 5 x 10<sup>-8</sup>) associated with folate levels...</p>

      <!-- Figure reference -->
      <fig id="fig1" position="float">
        <label>Figure 1</label>
        <caption>
          <title>Manhattan plot of GWAS results</title>
          <p>Genome-wide association results for serum folate levels.
            The red dashed line indicates genome-wide significance
            (P = 5 x 10<sup>-8</sup>).</p>
        </caption>
        <graphic xlink:href="ng.1234-F1.jpg" position="float"/>
      </fig>

      <!-- Table -->
      <table-wrap id="tbl1" position="float">
        <label>Table 1</label>
        <caption>
          <title>Top associated variants</title>
        </caption>
        <table frame="hsides" rules="groups">
          <thead>
            <tr>
              <th>rsID</th>
              <th>Gene</th>
              <th>P-value</th>
              <th>OR (95% CI)</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td>rs1801133</td>
              <td>MTHFR</td>
              <td>2.3 x 10<sup>-45</sup></td>
              <td>1.42 (1.35-1.49)</td>
            </tr>
            <tr>
              <td>rs1801131</td>
              <td>MTHFR</td>
              <td>8.7 x 10<sup>-23</sup></td>
              <td>1.18 (1.12-1.24)</td>
            </tr>
          </tbody>
        </table>
      </table-wrap>
    </sec>

    <sec id="sec4">
      <title>Discussion</title>
      <p>Our findings expand the genetic architecture of folate metabolism...</p>
    </sec>

  </body>

  <back>
    <!-- Supplementary materials -->
    <app-group>
      <app id="app1">
        <title>Supplementary Methods</title>
        <p>Detailed genotyping quality control procedures...</p>
      </app>
    </app-group>

    <!-- References -->
    <ref-list>
      <title>References</title>

      <ref id="ref1">
        <label>1</label>
        <element-citation publication-type="journal">
          <person-group person-group-type="author">
            <name>
              <surname>Frosst</surname>
              <given-names>P</given-names>
            </name>
            <name>
              <surname>Blom</surname>
              <given-names>HJ</given-names>
            </name>
          </person-group>
          <article-title>A candidate genetic risk factor for vascular
            disease: a common mutation in methylenetetrahydrofolate
            reductase</article-title>
          <source>Nat Genet</source>
          <year>1995</year>
          <volume>10</volume>
          <fpage>111</fpage>
          <lpage>113</lpage>
          <pub-id pub-id-type="pmid">7647779</pub-id>
          <pub-id pub-id-type="doi">10.1038/ng0595-111</pub-id>
        </element-citation>
      </ref>

      <ref id="ref15">
        <label>15</label>
        <element-citation publication-type="journal">
          <person-group person-group-type="author">
            <name>
              <surname>Chang</surname>
              <given-names>CC</given-names>
            </name>
          </person-group>
          <article-title>Second-generation PLINK</article-title>
          <source>GigaScience</source>
          <year>2015</year>
          <volume>4</volume>
          <fpage>7</fpage>
          <pub-id pub-id-type="pmid">25722852</pub-id>
        </element-citation>
      </ref>
    </ref-list>

  </back>
</article>
```

### 2.2 Section Structure Extraction

```typescript
interface JATSSection {
  id: string;
  title: string;
  type: 'introduction' | 'methods' | 'results' | 'discussion' |
        'conclusion' | 'supplementary' | 'other';
  content: string;
  subsections: JATSSection[];
  figures: JATSFigure[];
  tables: JATSTable[];
  references: string[];  // Reference IDs
}

// Section type detection
const SECTION_TYPE_PATTERNS: Record<string, RegExp> = {
  introduction: /^(introduction|background|overview)/i,
  methods: /^(methods?|materials?\s*(and|&)\s*methods?|experimental|procedures?)/i,
  results: /^(results?|findings?|outcomes?)/i,
  discussion: /^(discussion|interpretation)/i,
  conclusion: /^(conclusions?|summary|final\s*remarks?)/i,
};

function detectSectionType(title: string): string {
  for (const [type, pattern] of Object.entries(SECTION_TYPE_PATTERNS)) {
    if (pattern.test(title)) return type;
  }
  return 'other';
}
```

### 2.3 Figure/Table Handling

```typescript
interface JATSFigure {
  id: string;
  label: string;
  caption: string;
  title?: string;
  graphicUrl?: string;
  altText?: string;
}

interface JATSTable {
  id: string;
  label: string;
  caption: string;
  headers: string[];
  rows: string[][];
  footnotes?: string[];
}

// Extract table data for potential embedding
function tableToText(table: JATSTable): string {
  const lines: string[] = [];
  lines.push(`Table: ${table.label}. ${table.caption}`);
  lines.push(table.headers.join(' | '));
  for (const row of table.rows) {
    lines.push(row.join(' | '));
  }
  return lines.join('\n');
}
```

### 2.4 Citation/Reference Extraction

```typescript
interface JATSReference {
  id: string;
  label: string;
  authors: string[];
  title: string;
  source: string;  // Journal name
  year: string;
  volume?: string;
  pages?: string;
  pmid?: string;
  doi?: string;
  pmcid?: string;
}

// Build citation graph edges
function extractCitationEdges(
  articlePmid: string,
  references: JATSReference[]
): Array<{from: string; to: string; context?: string}> {
  return references
    .filter(ref => ref.pmid)
    .map(ref => ({
      from: articlePmid,
      to: ref.pmid!,
      context: ref.title
    }));
}
```

---

## 3. Metadata Standards

### 3.1 DOI System

```typescript
interface DOI {
  prefix: string;    // "10.1038" (publisher)
  suffix: string;    // "ng.1234" (article)
  full: string;      // "10.1038/ng.1234"
  url: string;       // "https://doi.org/10.1038/ng.1234"
}

// DOI resolution
const DOI_REGEX = /10\.\d{4,9}\/[-._;()/:A-Z0-9]+/i;

function parseDOI(input: string): DOI | null {
  const match = input.match(DOI_REGEX);
  if (!match) return null;

  const full = match[0];
  const [prefix, ...suffixParts] = full.split('/');

  return {
    prefix,
    suffix: suffixParts.join('/'),
    full,
    url: `https://doi.org/${full}`
  };
}
```

### 3.2 PMID vs PMCID vs DOI Relationships

```typescript
interface ArticleIdentifiers {
  pmid?: string;      // PubMed ID (e.g., "12345678")
  pmcid?: string;     // PubMed Central ID (e.g., "PMC12345678")
  doi?: string;       // DOI (e.g., "10.1038/ng.1234")
  pii?: string;       // Publisher Item ID

  // Relationships:
  // - PMID: Always unique, assigned by NLM
  // - PMCID: Only for PMC articles (open access subset)
  // - DOI: Publisher-assigned, may not exist for older articles
  // - One article can have all three
}

// ID conversion service (NCBI ID Converter)
// https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/
interface IDConverterResult {
  pmid: string;
  pmcid: string;
  doi: string;
  status: 'ok' | 'error';
  errmsg?: string;
}
```

### 3.3 MeSH Terms for Genetics Filtering

```typescript
// Key MeSH descriptors for genetics research
const GENETICS_MESH_TERMS = {
  // Core genetics
  polymorphism: {
    ui: 'D011110',
    name: 'Polymorphism, Genetic',
    children: [
      { ui: 'D020641', name: 'Polymorphism, Single Nucleotide' }
    ]
  },
  gwas: {
    ui: 'D055106',
    name: 'Genome-Wide Association Study'
  },
  pharmacogenetics: {
    ui: 'D010597',
    name: 'Pharmacogenetics'
  },

  // Gene-specific (examples)
  mthfr: {
    ui: 'D008780',
    name: 'Methylenetetrahydrofolate Reductase (NADPH2)'
  },

  // Qualifiers for filtering
  qualifiers: {
    genetics: { ui: 'Q000235', abbrev: 'GE' },
    metabolism: { ui: 'Q000378', abbrev: 'ME' },
    drug_effects: { ui: 'Q000494', abbrev: 'PD' }
  }
};

// Build PubMed query for genetics articles
function buildGeneticsQuery(geneSymbol: string): string {
  return `(${geneSymbol}[Title/Abstract] OR ${geneSymbol}[MeSH Terms])
    AND (polymorphism[MeSH Terms] OR variant[Title/Abstract]
    OR SNP[Title/Abstract] OR genetic[Title/Abstract])`;
}
```

### 3.4 Author Affiliation Data

```typescript
interface AuthorAffiliation {
  institution: string;
  department?: string;
  address?: {
    city?: string;
    state?: string;
    country: string;
    postal_code?: string;
  };
  email?: string;
  ror_id?: string;  // Research Organization Registry ID
}

// Affiliation parsing is challenging due to inconsistent formatting
// ROR (Research Organization Registry) provides standardized IDs
// https://ror.org/

function normalizeAffiliation(raw: string): AuthorAffiliation {
  // Basic parsing - production would use ML/NER
  const parts = raw.split(',').map(p => p.trim());

  return {
    institution: parts[0] || raw,
    address: {
      country: parts[parts.length - 1] || 'Unknown'
    }
  };
}
```

---

## 4. Text Extraction from PDFs

### 4.1 Tool Comparison

| Tool | Type | Quality | Speed | Features |
|------|------|---------|-------|----------|
| **PyMuPDF (fitz)** | Native | Good | Fast | Text + images, no ML |
| **pdfplumber** | Native | Good | Medium | Tables, precise layout |
| **GROBID** | ML-based | Excellent | Slow | Full structure, references |
| **Adobe API** | Cloud | Excellent | Variable | Best for scanned PDFs |
| **pdf.js** | JavaScript | Medium | Fast | Browser-compatible |

### 4.2 PyMuPDF Extraction

```python
import fitz  # PyMuPDF

def extract_pdf_pymupdf(pdf_path: str) -> dict:
    """Extract text and metadata from PDF using PyMuPDF."""
    doc = fitz.open(pdf_path)

    result = {
        'metadata': doc.metadata,
        'pages': [],
        'text': '',
        'images': []
    }

    for page_num, page in enumerate(doc):
        # Extract text with layout preservation
        text = page.get_text("text")

        # Extract text blocks with position info
        blocks = page.get_text("dict")["blocks"]

        # Extract images
        images = page.get_images()

        result['pages'].append({
            'number': page_num + 1,
            'text': text,
            'blocks': blocks,
            'image_count': len(images)
        })

        result['text'] += text + '\n\n'

    doc.close()
    return result
```

### 4.3 GROBID for Structured Extraction

```python
import requests

GROBID_URL = "http://localhost:8070/api/processFulltextDocument"

def extract_pdf_grobid(pdf_path: str) -> dict:
    """Extract structured content using GROBID ML service."""

    with open(pdf_path, 'rb') as f:
        response = requests.post(
            GROBID_URL,
            files={'input': f},
            data={
                'consolidateHeader': '1',
                'consolidateCitations': '1',
                'includeRawAffiliations': '1'
            }
        )

    # GROBID returns TEI XML
    tei_xml = response.text

    # Parse TEI XML to structured format
    return parse_grobid_tei(tei_xml)

def parse_grobid_tei(tei_xml: str) -> dict:
    """Parse GROBID TEI output to structured format."""
    from lxml import etree

    ns = {'tei': 'http://www.tei-c.org/ns/1.0'}
    root = etree.fromstring(tei_xml.encode())

    return {
        'title': root.find('.//tei:titleStmt/tei:title', ns).text,
        'abstract': root.find('.//tei:abstract', ns).text,
        'authors': extract_grobid_authors(root, ns),
        'sections': extract_grobid_sections(root, ns),
        'references': extract_grobid_references(root, ns)
    }
```

### 4.4 Quality Considerations

```typescript
interface PDFExtractionQuality {
  // Quality indicators
  characterAccuracy: number;       // OCR quality (0-1)
  structureConfidence: number;     // Section detection confidence
  tableExtraction: 'good' | 'partial' | 'failed';
  formulaExtraction: 'good' | 'partial' | 'failed';

  // Common issues
  issues: Array<
    | 'multi_column_merged'      // Two-column text merged incorrectly
    | 'headers_in_text'          // Page headers mixed with content
    | 'footnotes_misplaced'      // Footnotes in wrong location
    | 'tables_as_text'           // Tables converted to plain text
    | 'special_chars_garbled'    // Greek letters, math symbols wrong
    | 'citations_not_linked'     // Reference numbers not parsed
  >;
}

// Quality thresholds for acceptance
const QUALITY_THRESHOLDS = {
  minCharacterAccuracy: 0.95,
  minStructureConfidence: 0.8,
  maxIssueCount: 2
};

function assessExtractionQuality(extracted: any): PDFExtractionQuality {
  // Heuristic quality assessment
  const issues: string[] = [];

  // Check for multi-column merge issues
  if (extracted.text.match(/\w{50,}/)) {
    issues.push('multi_column_merged');
  }

  // Check for repeated headers
  const lines = extracted.text.split('\n');
  const lineFreq = new Map<string, number>();
  for (const line of lines) {
    const trimmed = line.trim();
    if (trimmed.length > 10 && trimmed.length < 100) {
      lineFreq.set(trimmed, (lineFreq.get(trimmed) || 0) + 1);
    }
  }
  for (const [, count] of lineFreq) {
    if (count > 3) {
      issues.push('headers_in_text');
      break;
    }
  }

  return {
    characterAccuracy: 0.98,  // Would use OCR confidence
    structureConfidence: 0.85,
    tableExtraction: 'partial',
    formulaExtraction: 'partial',
    issues: issues as any
  };
}
```

### 4.5 Structure Preservation Challenges

```typescript
// Challenge: Scientific papers have complex layouts
const LAYOUT_CHALLENGES = {
  // Two-column layout
  multiColumn: {
    problem: 'Text flows across columns incorrectly',
    solution: 'Use block coordinates to reorder text',
    tools: ['pdfplumber', 'GROBID']
  },

  // Figures and captions
  figureCaption: {
    problem: 'Captions separated from figures',
    solution: 'Match by proximity and numbering',
    tools: ['GROBID']
  },

  // Tables
  tables: {
    problem: 'Table structure lost in extraction',
    solution: 'Use table-specific extraction',
    tools: ['pdfplumber', 'Camelot', 'Tabula']
  },

  // Equations
  equations: {
    problem: 'Math rendered as images or garbled',
    solution: 'OCR with math recognition (Mathpix)',
    tools: ['Mathpix', 'LaTeX-OCR']
  },

  // References
  references: {
    problem: 'Citations not linked to bibliography',
    solution: 'Use GROBID or regex matching',
    tools: ['GROBID', 'AnyStyle']
  }
};
```

---

## 5. Embedding Considerations for RAG

### 5.1 Optimal Chunk Sizes

| Chunk Size | Tokens | Use Case | Pros | Cons |
|------------|--------|----------|------|------|
| **Sentence** | 20-50 | Precise retrieval | High precision | Lacks context |
| **Paragraph** | 100-200 | Balanced | Good context | May split ideas |
| **512 tokens** | 512 | Standard RAG | Efficient embedding | Fixed boundary |
| **1024 tokens** | 1024 | Long context | Full sections | Diluted relevance |
| **Section** | Variable | Semantic chunks | Natural boundaries | Size variance |

**Recommendation for Genetics Papers:**

```typescript
const CHUNK_CONFIG = {
  // Primary strategy: Section-based chunking with size limits
  strategy: 'hybrid',

  // For abstracts: Keep whole (usually < 512 tokens)
  abstract: {
    method: 'whole',
    maxTokens: 512
  },

  // For full text: Section-based with overlap
  fullText: {
    method: 'section',
    targetTokens: 512,
    maxTokens: 1024,
    overlapTokens: 50,  // Prevent context loss at boundaries
    splitOnParagraph: true
  },

  // For tables: Convert to text, embed separately
  tables: {
    method: 'separate',
    maxTokens: 256,
    includeCaption: true
  }
};
```

### 5.2 What to Embed

```typescript
interface EmbeddableContent {
  // Always embed
  abstract: string;

  // Embed if available
  sections?: {
    introduction: string;
    methods: string;
    results: string;
    discussion: string;
    conclusion: string;
  };

  // Selective embedding
  figures?: Array<{
    caption: string;
    description?: string;
  }>;

  tables?: Array<{
    caption: string;
    content: string;  // Converted to text
  }>;

  // Skip embedding (metadata only)
  references: never;  // Use graph links instead
  authorInfo: never;  // Store as metadata
}

// Embedding priority for storage optimization
enum EmbeddingPriority {
  REQUIRED = 1,    // Abstract (always embed)
  HIGH = 2,        // Results, conclusions (embed if space)
  MEDIUM = 3,      // Introduction, discussion
  LOW = 4,         // Methods, supplementary
  SKIP = 5         // References, acknowledgments
}
```

### 5.3 Embedding Strategies by Article Type

```typescript
const EMBEDDING_BY_TYPE: Record<string, EmbeddingStrategy> = {
  // Research articles: Full content
  'research-article': {
    sections: ['abstract', 'introduction', 'results', 'discussion', 'conclusion'],
    chunkSize: 512,
    embedTables: true,
    embedFigureCaptions: true
  },

  // Reviews: Dense content, larger chunks
  'review-article': {
    sections: ['abstract', 'all'],  // Embed all sections
    chunkSize: 1024,  // Larger for comprehensive context
    embedTables: true,
    embedFigureCaptions: true
  },

  // Case reports: Short, full embed
  'case-report': {
    sections: ['abstract', 'case', 'discussion'],
    chunkSize: 512,
    embedTables: false,
    embedFigureCaptions: true
  },

  // Meta-analyses: Focus on results
  'meta-analysis': {
    sections: ['abstract', 'results', 'forest-plots'],
    chunkSize: 512,
    embedTables: true,  // Summary tables important
    embedFigureCaptions: true
  }
};
```

### 5.4 Metadata to Preserve Alongside Embeddings

```typescript
interface PaperChunkMetadata {
  // Required identifiers
  pmid: string;
  chunkId: string;           // "{pmid}-{section}-{chunk_num}"

  // Source location
  section: string;           // "abstract", "results", etc.
  chunkIndex: number;        // Position within section
  pageNumbers?: number[];    // Original PDF pages

  // Content type
  contentType: 'text' | 'table' | 'figure_caption';

  // Article metadata (denormalized for filtering)
  publicationYear: number;
  journal: string;
  articleType: string;

  // Genetics-specific
  mentionedGenes: string[];  // Gene symbols found in chunk
  mentionedSNPs: string[];   // RS numbers found in chunk
  meshTerms: string[];       // MeSH descriptors

  // Quality indicators
  embeddingModel: string;    // "all-MiniLM-L6-v2"
  tokenCount: number;

  // For citations
  relatedPMIDs?: string[];   // Papers this chunk cites
}
```

### 5.5 RuVector Schema for Papers

```typescript
import { RuVector } from 'ruvector';

// Paper collection with chunks
const paperSchema = {
  collections: {
    // Main article metadata
    articles: {
      dimension: 384,
      distanceMetric: 'cosine',
      compression: 'auto',
      properties: {
        pmid: { type: 'number', indexed: true },
        pmcid: { type: 'string', indexed: true },
        doi: { type: 'string', indexed: true },
        title: { type: 'string' },
        journal: { type: 'string', indexed: true },
        publication_year: { type: 'number', indexed: true },
        article_type: { type: 'string', indexed: true },
        mesh_terms: { type: 'array' },
        authors: { type: 'array' },
        cited_genes: { type: 'array', indexed: true },
        cited_snps: { type: 'array', indexed: true }
      }
    },

    // Individual chunks for RAG
    paper_chunks: {
      dimension: 384,
      distanceMetric: 'cosine',
      compression: 'auto',
      properties: {
        article_id: { type: 'string', indexed: true },  // Link to articles
        pmid: { type: 'number', indexed: true },
        section: { type: 'string', indexed: true },
        chunk_index: { type: 'number' },
        content_type: { type: 'string', indexed: true },
        token_count: { type: 'number' },
        mentioned_genes: { type: 'array', indexed: true },
        mentioned_snps: { type: 'array', indexed: true }
      }
    }
  },

  graph: {
    enabled: true,
    relationships: [
      // Article-level relationships
      'CITES',           // (:Article)-[:CITES]->(:Article)
      'MENTIONS_GENE',   // (:Article)-[:MENTIONS_GENE]->(:Gene)
      'MENTIONS_SNP',    // (:Article)-[:MENTIONS_SNP]->(:SNP)
      'PUBLISHED_IN',    // (:Article)-[:PUBLISHED_IN]->(:Journal)

      // Chunk relationships
      'CHUNK_OF',        // (:PaperChunk)-[:CHUNK_OF]->(:Article)
      'NEXT_CHUNK',      // (:PaperChunk)-[:NEXT_CHUNK]->(:PaperChunk)

      // Cross-entity
      'CITED_IN'         // (:SNP)-[:CITED_IN]->(:Article)
    ]
  }
};
```

---

## 6. Example Data Structures (TypeScript)

### 6.1 Core Paper Interfaces

```typescript
// Full paper representation
interface ResearchPaper {
  // Identifiers
  pmid: number;
  pmcid?: string;
  doi?: string;
  pii?: string;

  // Basic metadata
  title: string;
  abstract: string;
  journal: Journal;
  publicationDate: Date;
  articleType: ArticleType;
  language: string;

  // Authors
  authors: Author[];

  // Content (PMC only)
  fullText?: FullTextContent;

  // Controlled vocabulary
  meshTerms: MeSHTerm[];
  keywords: string[];
  chemicals: Chemical[];

  // Extracted entities
  mentionedGenes: string[];
  mentionedSNPs: string[];
  mentionedDiseases: string[];

  // Embeddings (computed)
  abstractEmbedding?: Float32Array;
  chunks?: PaperChunk[];

  // Graph relationships
  citations: string[];       // PMIDs this paper cites
  citedBy?: string[];        // PMIDs that cite this paper
}

interface Journal {
  title: string;
  isoAbbreviation: string;
  issn?: string;
  eIssn?: string;
  volume?: string;
  issue?: string;
  pages?: string;
}

interface Author {
  lastName: string;
  foreName?: string;
  initials: string;
  orcid?: string;
  affiliations: Affiliation[];
  isCorresponding?: boolean;
  equalContribution?: boolean;
}

interface Affiliation {
  institution: string;
  department?: string;
  city?: string;
  country: string;
  email?: string;
  rorId?: string;
}

type ArticleType =
  | 'journal-article'
  | 'review'
  | 'meta-analysis'
  | 'case-report'
  | 'clinical-trial'
  | 'letter'
  | 'editorial'
  | 'preprint';

interface MeSHTerm {
  descriptorUI: string;
  descriptorName: string;
  qualifierUI?: string;
  qualifierName?: string;
  majorTopic: boolean;
}

interface Chemical {
  registryNumber: string;
  name: string;
  ui?: string;
}
```

### 6.2 Full-Text Content Structure

```typescript
interface FullTextContent {
  sections: Section[];
  figures: Figure[];
  tables: Table[];
  equations: Equation[];
  supplementaryMaterials?: SupplementaryMaterial[];
  references: Reference[];
}

interface Section {
  id: string;
  title: string;
  type: SectionType;
  content: string;
  subsections: Section[];

  // Positions of inline elements
  figureRefs: { id: string; position: number }[];
  tableRefs: { id: string; position: number }[];
  citationRefs: { id: string; position: number }[];
}

type SectionType =
  | 'introduction'
  | 'methods'
  | 'results'
  | 'discussion'
  | 'conclusion'
  | 'abstract'
  | 'supplementary'
  | 'other';

interface Figure {
  id: string;
  label: string;
  caption: string;
  graphicUrl?: string;
  altText?: string;
  panels?: FigurePanel[];
}

interface FigurePanel {
  label: string;  // "A", "B", etc.
  description: string;
}

interface Table {
  id: string;
  label: string;
  caption: string;
  headers: TableHeader[];
  rows: TableRow[];
  footnotes?: string[];
}

interface TableHeader {
  text: string;
  colspan?: number;
  rowspan?: number;
}

interface TableRow {
  cells: TableCell[];
  isHeader?: boolean;
}

interface TableCell {
  text: string;
  colspan?: number;
  rowspan?: number;
}

interface Equation {
  id: string;
  latex?: string;
  mathml?: string;
  plainText?: string;
}

interface Reference {
  id: string;
  label: string;
  authors: string[];
  title: string;
  source: string;
  year: number;
  volume?: string;
  pages?: string;
  pmid?: number;
  doi?: string;
  pmcid?: string;
}
```

### 6.3 Chunk Representation

```typescript
interface PaperChunk {
  id: string;                    // Unique chunk ID
  articlePmid: number;           // Parent article

  // Location
  section: SectionType;
  chunkIndex: number;

  // Content
  text: string;
  tokenCount: number;
  contentType: 'text' | 'table' | 'figure_caption';

  // Embedding
  embedding: Float32Array;
  embeddingModel: string;

  // Extracted entities (for graph links)
  mentionedGenes: string[];
  mentionedSNPs: string[];

  // For reconstruction
  previousChunkId?: string;
  nextChunkId?: string;
}
```

### 6.4 RuVector Storage Models

```typescript
// RuVector article node
interface ArticleNode {
  id: string;                    // "pmid:12345678"
  pmid: number;
  pmcid?: string;
  doi?: string;

  title: string;
  abstract: string;
  journal: string;
  publication_year: number;
  article_type: ArticleType;

  // For filtering
  mesh_terms: string[];          // MeSH descriptor names
  mentioned_genes: string[];     // Gene symbols
  mentioned_snps: string[];      // RS numbers

  // Embedding of title + abstract
  embedding: Float32Array;
}

// RuVector chunk node
interface ChunkNode {
  id: string;                    // "pmid:12345678:results:0"
  article_id: string;            // Link to parent
  pmid: number;

  section: string;
  chunk_index: number;
  content_type: string;

  text: string;
  token_count: number;

  mentioned_genes: string[];
  mentioned_snps: string[];

  embedding: Float32Array;
}
```

### 6.5 Graph Relationships

```cypher
// Article citations
(:Article {pmid: 12345678})-[:CITES {context: "Methods section"}]->
  (:Article {pmid: 98765432})

// Article mentions entities
(:Article {pmid: 12345678})-[:MENTIONS_GENE {
  section: "results",
  count: 15,
  context: "MTHFR C677T polymorphism was associated..."
}]->(:Gene {symbol: "MTHFR"})

(:Article {pmid: 12345678})-[:MENTIONS_SNP {
  section: "results",
  count: 23,
  context: "rs1801133 (C677T)..."
}]->(:SNP {rs_number: "rs1801133"})

// SNP evidence from papers
(:SNP {rs_number: "rs1801133"})-[:CITED_IN {
  evidence_level: "established",
  study_type: "meta-analysis",
  sample_size: 50000
}]->(:Article {pmid: 12345678})

// Chunk navigation
(:PaperChunk {id: "pmid:12345678:results:0"})-[:NEXT_CHUNK]->
  (:PaperChunk {id: "pmid:12345678:results:1"})

(:PaperChunk)-[:CHUNK_OF]->(:Article)
```

### 6.6 Query Examples

```typescript
// Find papers about a specific SNP
async function findPapersForSNP(rsNumber: string): Promise<ArticleNode[]> {
  const result = await db.cypher(`
    MATCH (snp:SNP {rs_number: $rsNumber})<-[:MENTIONS_SNP]-(article:Article)
    RETURN article
    ORDER BY article.publication_year DESC
    LIMIT 50
  `, { rsNumber });

  return result.records.map(r => r.get('article'));
}

// RAG: Find relevant chunks for a question
async function findRelevantChunks(
  question: string,
  filters?: { genes?: string[]; snps?: string[] }
): Promise<ChunkNode[]> {
  const embedding = await db.embed(question);

  let query = `
    MATCH (chunk:PaperChunk)
    WHERE chunk.embedding <-> $embedding < 0.4
  `;

  if (filters?.genes?.length) {
    query += `
      AND ANY(gene IN chunk.mentioned_genes WHERE gene IN $genes)
    `;
  }

  if (filters?.snps?.length) {
    query += `
      AND ANY(snp IN chunk.mentioned_snps WHERE snp IN $snps)
    `;
  }

  query += `
    RETURN chunk,
           (1 - chunk.embedding <-> $embedding) as similarity
    ORDER BY similarity DESC
    LIMIT 10
  `;

  const result = await db.cypher(query, {
    embedding,
    genes: filters?.genes || [],
    snps: filters?.snps || []
  });

  return result.records.map(r => ({
    ...r.get('chunk'),
    similarity: r.get('similarity')
  }));
}

// Build citation graph for a paper
async function getCitationNetwork(pmid: number, depth: number = 2): Promise<any> {
  const result = await db.cypher(`
    MATCH path = (start:Article {pmid: $pmid})-[:CITES*1..${depth}]->(cited:Article)
    RETURN path
  `, { pmid });

  return result;
}

// Find papers connecting SNP to phenotype
async function findEvidencePath(
  rsNumber: string,
  phenotypeName: string
): Promise<any> {
  const result = await db.cypher(`
    MATCH (snp:SNP {rs_number: $rsNumber}),
          (phenotype:Phenotype {name: $phenotypeName})
    MATCH path = (snp)<-[:MENTIONS_SNP]-(article:Article)
                 -[:MENTIONS_PHENOTYPE]->(phenotype)
    RETURN path, article
    ORDER BY article.publication_year DESC
  `, { rsNumber, phenotypeName });

  return result;
}
```

---

## 7. Data Pipeline Architecture

### 7.1 Ingestion Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                     PAPER INGESTION PIPELINE                         │
└─────────────────────────────────────────────────────────────────────┘
                                  │
           ┌──────────────────────┼──────────────────────┐
           │                      │                      │
           ▼                      ▼                      ▼
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│  PubMed Baseline │    │  PubMed Daily   │    │  PMC Open       │
│  (Annual XML)    │    │  Updates        │    │  Access         │
│  ~36M articles   │    │  ~1K/day        │    │  ~8M articles   │
└─────────────────┘    └─────────────────┘    └─────────────────┘
           │                      │                      │
           └──────────────────────┼──────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │     XML Parser          │
                    │  (Streaming for large   │
                    │   files)                │
                    └─────────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │   Genetics Filter       │
                    │  - MeSH: Polymorphism,  │
                    │    Pharmacogenetics     │
                    │  - Keywords: SNP, gene, │
                    │    variant, GWAS        │
                    │  - Reduces to ~2M       │
                    └─────────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │   Entity Extraction     │
                    │  - Gene symbols (NER)   │
                    │  - RS numbers (regex)   │
                    │  - Disease names (NER)  │
                    └─────────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │    Chunking Engine      │
                    │  - Abstract: whole      │
                    │  - Full text: sections  │
                    │  - Tables: converted    │
                    │  - Target: 512 tokens   │
                    └─────────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────┐
                    │   Embedding Generation  │
                    │  - Model: MiniLM-L6-v2  │
                    │  - Batch size: 100      │
                    │  - GPU: 15K/sec         │
                    └─────────────────────────┘
                                  │
                    ┌─────────────┴─────────────┐
                    │                           │
                    ▼                           ▼
          ┌─────────────────┐        ┌─────────────────┐
          │    RuVector     │        │   PostgreSQL    │
          │  - Articles     │        │  - Full text    │
          │  - Chunks       │        │  - User reports │
          │  - Graph edges  │        │  - Audit logs   │
          └─────────────────┘        └─────────────────┘
```

### 7.2 Incremental Update Strategy

```typescript
interface UpdateStrategy {
  // Daily PubMed updates
  pubmedDaily: {
    source: 'ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/',
    schedule: '0 4 * * *',  // 4 AM daily
    process: 'incremental'
  };

  // Weekly PMC updates
  pmcWeekly: {
    source: 'ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/',
    schedule: '0 2 * * 0',  // 2 AM Sunday
    process: 'incremental'
  };

  // Entity re-extraction (monthly)
  entityRefresh: {
    schedule: '0 0 1 * *',  // 1st of month
    process: 'full'
  };

  // Embedding regeneration (on model update)
  embeddingRefresh: {
    trigger: 'model_version_change',
    process: 'full'
  };
}
```

---

## 8. Size Estimates for Papers

### 8.1 Storage Requirements

| Component | Count | Size per Item | Total |
|-----------|-------|---------------|-------|
| **Articles (metadata + abstract)** | 2M | ~5 KB | ~10 GB |
| **Article embeddings (384d)** | 2M | 1.5 KB (f32) | ~3 GB |
| **Article embeddings (RuVector compressed)** | 2M | ~400 bytes | ~800 MB |
| **Chunks** | 10M | ~1 KB | ~10 GB |
| **Chunk embeddings (384d)** | 10M | 1.5 KB (f32) | ~15 GB |
| **Chunk embeddings (RuVector compressed)** | 10M | ~400 bytes | ~4 GB |
| **Full text (PMC subset)** | 500K | ~50 KB | ~25 GB |
| **Graph edges** | 50M | ~100 bytes | ~5 GB |

### 8.2 Estimated Totals

| Configuration | Articles | Chunks | RuVector Storage |
|---------------|----------|--------|------------------|
| **MVP** | 500K | 2M | ~2 GB |
| **Standard** | 2M | 10M | ~6 GB |
| **Comprehensive** | 5M | 25M | ~15 GB |

---

## 9. Implementation Checklist

### Phase 1: PubMed Abstracts
- [ ] Set up PubMed baseline download pipeline
- [ ] Implement streaming XML parser
- [ ] Create genetics relevance filter (MeSH + keywords)
- [ ] Extract gene/SNP mentions with NER
- [ ] Generate abstract embeddings
- [ ] Build citation graph edges
- [ ] Load into RuVector articles collection

### Phase 2: Full-Text (PMC)
- [ ] Download PMC Open Access subset
- [ ] Implement JATS XML parser
- [ ] Create section-based chunking
- [ ] Extract tables and figure captions
- [ ] Generate chunk embeddings
- [ ] Build chunk navigation graph
- [ ] Link chunks to entities

### Phase 3: PDF Fallback
- [ ] Set up GROBID service
- [ ] Implement PDF extraction pipeline
- [ ] Quality assessment and filtering
- [ ] Integration with main pipeline

### Phase 4: RAG Integration
- [ ] Implement hybrid search (vector + graph)
- [ ] Build context construction for Claude
- [ ] Create citation tracking for answers
- [ ] Add relevance feedback loop

---

## References

- PubMed DTD: https://dtd.nlm.nih.gov/ncbi/pubmed/out/
- JATS XML: https://jats.nlm.nih.gov/
- MeSH Browser: https://meshb.nlm.nih.gov/
- GROBID: https://github.com/kermitt2/grobid
- PyMuPDF: https://pymupdf.readthedocs.io/
- RuVector: https://ruvector.io/

---

*This document is part of the Gene Knowledge Base technical documentation.*
