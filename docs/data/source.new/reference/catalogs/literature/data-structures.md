---
id: literature-data-structures
title: Literature Data Structures for Genetics Knowledge Base
category: literature
tier: 2
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [literature, data-structures, xml, jats, ruvector, schema]
---

# Research Paper Data Structures for Genetics Knowledge Base

**Last Updated:** January 2026
**Purpose:** Data formats, parsing strategies, and embedding approaches for biomedical research papers
**Target System:** RuVector (vector database with Cypher graph queries)
**Parent:** [../_index.md](../_index.md)

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

PubMed uses the `pubmed_230101.dtd` Document Type Definition (updated annually).

```xml
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMed 230101//EN"
    "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_230101.dtd">
```

### 1.2 Key XML Elements

```xml
<PubmedArticleSet>
  <PubmedArticle>
    <MedlineCitation>
      <PMID Version="1">12345678</PMID>
      <Article>
        <Journal>
          <Title>Nature Genetics</Title>
          <ISSN>1061-4036</ISSN>
        </Journal>
        <ArticleTitle>MTHFR C677T polymorphism analysis</ArticleTitle>
        <Abstract>
          <AbstractText Label="BACKGROUND">...</AbstractText>
          <AbstractText Label="METHODS">...</AbstractText>
          <AbstractText Label="RESULTS">...</AbstractText>
        </Abstract>
        <AuthorList>
          <Author>
            <LastName>Smith</LastName>
            <ForeName>John A</ForeName>
          </Author>
        </AuthorList>
      </Article>
      <MeshHeadingList>
        <MeshHeading>
          <DescriptorName>Polymorphism, Single Nucleotide</DescriptorName>
          <QualifierName>genetics</QualifierName>
        </MeshHeading>
      </MeshHeadingList>
      <KeywordList>
        <Keyword>SNP</Keyword>
        <Keyword>MTHFR</Keyword>
      </KeywordList>
    </MedlineCitation>
  </PubmedArticle>
</PubmedArticleSet>
```

### 1.3 Available Metadata Fields

- PMID (unique identifier)
- Title, Abstract (with structured labels)
- Authors (names, affiliations, ORCID)
- Journal info (name, ISSN, volume, issue, pages)
- Publication date
- MeSH terms (Medical Subject Headings)
- Keywords, Publication type
- DOI, PMC ID links, Language
- References, Grant information

---

## 2. PMC Full-Text Format (JATS XML)

### 2.1 JATS XML Schema

PMC uses the Journal Article Tag Suite (JATS) standard:

```xml
<article xmlns:xlink="http://www.w3.org/1999/xlink" article-type="research-article">
  <front>
    <article-meta>
      <article-id pub-id-type="pmid">12345678</article-id>
      <article-id pub-id-type="pmc">PMC1234567</article-id>
      <article-id pub-id-type="doi">10.1038/ng.1234</article-id>
      <title-group>
        <article-title>Genome-wide association study</article-title>
      </title-group>
      <abstract>
        <sec><title>Background</title><p>...</p></sec>
        <sec><title>Results</title><p>...</p></sec>
      </abstract>
    </article-meta>
  </front>
  <body>
    <sec id="sec1">
      <title>Introduction</title>
      <p>MTHFR plays a critical role...</p>
    </sec>
    <sec id="sec2">
      <title>Methods</title>
      <sec id="sec2-1">
        <title>Study Population</title>
        <p>We recruited 50,000 participants...</p>
      </sec>
    </sec>
    <sec id="sec3">
      <title>Results</title>
      <table-wrap id="tbl1">
        <label>Table 1</label>
        <caption><title>Top associated variants</title></caption>
        <table>
          <thead>
            <tr><th>rsID</th><th>Gene</th><th>P-value</th></tr>
          </thead>
          <tbody>
            <tr><td>rs1801133</td><td>MTHFR</td><td>2.3e-45</td></tr>
          </tbody>
        </table>
      </table-wrap>
    </sec>
  </body>
  <back>
    <ref-list>
      <ref id="ref1">
        <element-citation publication-type="journal">
          <person-group><name><surname>Frosst</surname></name></person-group>
          <article-title>A candidate genetic risk factor</article-title>
          <source>Nat Genet</source>
          <year>1995</year>
          <pub-id pub-id-type="pmid">7647779</pub-id>
        </element-citation>
      </ref>
    </ref-list>
  </back>
</article>
```

### 2.2 Full-Text Elements

- Complete article body with sections
- Tables, figures (with captions)
- Supplementary materials
- References with full citation data
- Author affiliations, Funding information
- Conflict of interest statements

---

## 3. Core Data Structures (TypeScript)

### 3.1 Research Paper Interface

```typescript
interface ResearchPaper {
  // Identifiers
  pmid: number;
  pmcid?: string;
  doi?: string;

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

interface Author {
  lastName: string;
  foreName?: string;
  initials: string;
  orcid?: string;
  affiliations: Affiliation[];
}

interface MeSHTerm {
  descriptorUI: string;
  descriptorName: string;
  qualifierName?: string;
  majorTopic: boolean;
}
```

### 3.2 Full-Text Content Structure

```typescript
interface FullTextContent {
  sections: Section[];
  figures: Figure[];
  tables: Table[];
  references: Reference[];
}

interface Section {
  id: string;
  title: string;
  type: 'introduction' | 'methods' | 'results' | 'discussion' | 'conclusion';
  content: string;
  subsections: Section[];
}

interface Table {
  id: string;
  label: string;
  caption: string;
  headers: string[];
  rows: string[][];
}

interface Reference {
  id: string;
  authors: string[];
  title: string;
  source: string;
  year: number;
  pmid?: number;
  doi?: string;
}
```

### 3.3 Chunk Representation

```typescript
interface PaperChunk {
  id: string;                    // "pmid:12345678:results:0"
  articlePmid: number;

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

  // Extracted entities
  mentionedGenes: string[];
  mentionedSNPs: string[];

  // For reconstruction
  previousChunkId?: string;
  nextChunkId?: string;
}
```

---

## 4. RuVector Storage Schema

### 4.1 Collections

```typescript
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
        mesh_terms: { type: 'array' },
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
        article_id: { type: 'string', indexed: true },
        pmid: { type: 'number', indexed: true },
        section: { type: 'string', indexed: true },
        chunk_index: { type: 'number' },
        mentioned_genes: { type: 'array', indexed: true },
        mentioned_snps: { type: 'array', indexed: true }
      }
    }
  }
};
```

### 4.2 Graph Relationships

```cypher
// Article citations
(:Article {pmid: 12345678})-[:CITES]->(:Article {pmid: 98765432})

// Article mentions entities
(:Article)-[:MENTIONS_GENE {section: "results"}]->(:Gene {symbol: "MTHFR"})
(:Article)-[:MENTIONS_SNP {count: 23}]->(:SNP {rs_number: "rs1801133"})

// SNP evidence from papers
(:SNP {rs_number: "rs1801133"})-[:CITED_IN {
  evidence_level: "established",
  study_type: "meta-analysis"
}]->(:Article)

// Chunk navigation
(:PaperChunk)-[:NEXT_CHUNK]->(:PaperChunk)
(:PaperChunk)-[:CHUNK_OF]->(:Article)
```

---

## 5. Embedding Considerations for RAG

### 5.1 Optimal Chunk Sizes

| Chunk Size | Tokens | Use Case | Pros | Cons |
|------------|--------|----------|------|------|
| **Sentence** | 20-50 | Precise retrieval | High precision | Lacks context |
| **Paragraph** | 100-200 | Balanced | Good context | May split ideas |
| **512 tokens** | 512 | Standard RAG | Efficient | Fixed boundary |
| **Section** | Variable | Semantic chunks | Natural boundaries | Size variance |

### 5.2 What to Embed

```typescript
interface EmbeddableContent {
  // Always embed
  abstract: string;

  // Embed if available
  sections?: {
    introduction: string;
    results: string;
    discussion: string;
    conclusion: string;
  };

  // Selective embedding
  figures?: Array<{ caption: string }>;
  tables?: Array<{ caption: string; content: string }>;
}

enum EmbeddingPriority {
  REQUIRED = 1,    // Abstract (always)
  HIGH = 2,        // Results, conclusions
  MEDIUM = 3,      // Introduction, discussion
  LOW = 4,         // Methods, supplementary
  SKIP = 5         // References, acknowledgments
}
```

### 5.3 Metadata Alongside Embeddings

```typescript
interface PaperChunkMetadata {
  // Identifiers
  pmid: string;
  chunkId: string;

  // Source location
  section: string;
  chunkIndex: number;
  pageNumbers?: number[];

  // Article metadata (denormalized)
  publicationYear: number;
  journal: string;
  articleType: string;

  // Genetics-specific
  mentionedGenes: string[];
  mentionedSNPs: string[];
  meshTerms: string[];

  // Quality indicators
  embeddingModel: string;
  tokenCount: number;

  // For citations
  relatedPMIDs?: string[];
}
```

---

## 6. Query Examples

### 6.1 Find Papers for SNP

```typescript
async function findPapersForSNP(rsNumber: string): Promise<ArticleNode[]> {
  const result = await db.cypher(`
    MATCH (snp:SNP {rs_number: $rsNumber})<-[:MENTIONS_SNP]-(article:Article)
    RETURN article
    ORDER BY article.publication_year DESC
    LIMIT 50
  `, { rsNumber });

  return result.records.map(r => r.get('article'));
}
```

### 6.2 RAG: Find Relevant Chunks

```typescript
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

  query += `
    RETURN chunk,
           (1 - chunk.embedding <-> $embedding) as similarity
    ORDER BY similarity DESC
    LIMIT 10
  `;

  return await db.cypher(query, {
    embedding,
    genes: filters?.genes || []
  });
}
```

### 6.3 Citation Network

```typescript
async function getCitationNetwork(pmid: number, depth: number = 2) {
  return await db.cypher(`
    MATCH path = (start:Article {pmid: $pmid})-[:CITES*1..${depth}]->(cited:Article)
    RETURN path
  `, { pmid });
}
```

---

## 7. Size Estimates

### 7.1 Storage per Paper

| Component | Abstract | Full-Text |
|-----------|----------|-----------|
| **Text only** | 1.5 KB | 30-50 KB |
| **+ Metadata** | 2 KB | 35-55 KB |
| **+ Embedding (384d)** | 3.5 KB | 36-56 KB |
| **+ 5 chunk embeddings** | 9 KB | 60-100 KB |

### 7.2 Database Size (5M Papers)

| Component | Size |
|-----------|------|
| **Articles (metadata + abstract)** | ~10 GB |
| **Article embeddings (compressed)** | ~800 MB |
| **Chunks** | ~10 GB |
| **Chunk embeddings (compressed)** | ~4 GB |
| **Graph edges** | ~5 GB |
| **Total** | **~30 GB** |

---

## 8. Implementation Priorities

### Phase 1: PubMed Abstracts
- Implement PubMed XML parser
- Extract title, abstract, metadata
- Generate embeddings
- Load into RuVector articles collection
- Build CITED_IN relationships

### Phase 2: Full-Text (PMC)
- Implement JATS XML parser
- Extract sections, tables, figures
- Create chunks
- Generate chunk embeddings
- Build chunk navigation graph

### Phase 3: RAG Integration
- Implement hybrid search (vector + graph)
- Build context construction
- Add citation tracking
- Relevance feedback loop

---

## Download

| Source | Method | URL/Command |
|--------|--------|-------------|
| **PubMed XML** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/` |
| **PMC JATS** | FTP | `ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/` |
| **Annotations** | Europe PMC | `https://europepmc.org/downloads/annotations` |

**Access Requirements:** All freely accessible; API key recommended for higher rate limits.

## Data Format

| Format | Description |
|--------|-------------|
| Primary | XML (PubMed DTD, JATS) |
| Alternative | JSON, Parquet |
| Embeddings | Float32 arrays |
| Graph | RuVector/Neo4j Cypher |
| Encoding | UTF-8 |

## Schema

### Core Fields

| Field | Type | Description | Example |
|-------|------|-------------|---------|
| `pmid` | integer | PubMed identifier | 12345678 |
| `sections` | object | Parsed JATS sections | {"intro": "...", "methods": "..."} |
| `chunks` | array | Text chunks for RAG | [{text, embedding, position}] |
| `entities` | array | Extracted NER entities | ["rs1801133", "MTHFR"] |

### Relationships

| Relation | Target | Cardinality |
|----------|--------|-------------|
| `has_section` | Section | 1:N |
| `cites` | Article | N:M |
| `mentions` | Entity | N:M |

## Sample Data

### Example Structured Article
```json
{
  "pmid": 12345678,
  "pmcid": "PMC1234567",
  "sections": {
    "title": "MTHFR polymorphisms...",
    "abstract": "Background: ...",
    "introduction": "The MTHFR gene...",
    "methods": "We conducted..."
  },
  "chunks": [
    {"text": "Background: The MTHFR...", "embedding": [0.023, ...], "section": "abstract"}
  ]
}
```

### Sample Query Result
| pmid | section | chunk_count | entities |
|------|---------|-------------|----------|
| 12345678 | methods | 15 | ["MTHFR", "rs1801133"] |
| 23456789 | results | 22 | ["CYP2D6", "codeine"] |

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| PubMed | Public domain | Yes |
| PMC OA | CC BY/CC0 | Yes (OA subset) |
| Europe PMC | CC BY | Yes |

## Data Set Size

| Metric | Value |
|--------|-------|
| PubMed records | 39M+ citations |
| PMC full-text | 3.4M+ articles |
| Chunks per article | ~20 average |
| Total chunks | ~70M estimated |
| Last updated | January 2026 |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `JATS XML` | Journal Article Tag Suite - NISO standard for full-text article markup | PMC article format |
| `PubMed DTD` | Document Type Definition specifying PubMed XML structure | pubmed_230101.dtd |
| `PMID` | PubMed Identifier - unique integer for each citation | PMID: 12345678 |
| `PMCID` | PubMed Central Identifier - unique ID for full-text articles | PMC1234567 |
| `embedding` | Dense vector representation for semantic similarity search | Float32Array of 384 values |
| `chunk` | Text segment sized for embedding and RAG retrieval | 512-token section |
| `Cypher` | Graph query language used with Neo4j and RuVector | MATCH (n)-[:CITES]->(m) |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| `RuVector` | Vector database with hybrid vector + graph capabilities | Storage backend |
| `collection` | A set of vectors with shared schema in RuVector | articles, paper_chunks |
| `cosine distance` | Similarity metric measuring angle between vectors | Distance metric |
| `MeSH terms` | Medical Subject Headings - controlled vocabulary for PubMed | Article indexing |
| `CITED_IN relationship` | Graph edge linking SNPs/genes to citing articles | Graph schema |
| `tiered compression` | Progressive lossy compression (f32 to PQ4) for storage efficiency | Storage optimization |
| `section-based chunking` | Splitting articles by semantic sections (intro, methods, results) | Chunking strategy |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| JATS | Journal Article Tag Suite | Full-text XML standard |
| DTD | Document Type Definition | XML schema format |
| PMID | PubMed Identifier | Citation unique ID |
| PMC | PubMed Central | Full-text archive |
| DOI | Digital Object Identifier | Persistent paper ID |
| MeSH | Medical Subject Headings | Controlled vocabulary |
| RAG | Retrieval-Augmented Generation | Search + LLM pattern |
| PQ | Product Quantization | Vector compression |
| HNSW | Hierarchical Navigable Small World | Vector index algorithm |
| NER | Named Entity Recognition | Entity extraction |
| ORCID | Open Researcher and Contributor ID | Author identifier |

---

## References

- [PubMed DTD](https://dtd.nlm.nih.gov/ncbi/pubmed/out/)
- [JATS XML](https://jats.nlm.nih.gov/)
- [MeSH Browser](https://meshb.nlm.nih.gov/)
- [RuVector](https://ruvector.io/)

---

*This document is part of the Gene Knowledge Base technical documentation.*

**Note:** The original data-structures.md file contains extensive TypeScript interface definitions and XML examples. This summary focuses on key structures. Refer to the original file in `/docs/data/source/literature/` for complete type definitions.
