---
id: reference-literature-formats
title: Literature Data Structures for Genetics Knowledge Base
category: reference
parent: _index.md
last_updated: 2026-01-23
status: active
migrated_from: databases/literature/data-structures.md
tags: [literature, data-structures, xml, jats, ruvector, schema]
---

# Research Paper Data Structures for Genetics Knowledge Base

**Last Updated:** January 2026
**Purpose:** Data formats, parsing strategies, and embedding approaches for biomedical research papers
**Target System:** RuVector (vector database with Cypher graph queries)
**Parent:** [Format Specifications](./_index.md)

---

## Overview

Research papers are critical for the genetics knowledge base, providing:
1. **Evidence** - Citations for SNP-phenotype associations
2. **Context** - Detailed mechanisms and study findings
3. **Discovery** - New relationships via semantic search and RAG
4. **Authority** - Academic credibility for recommendations

---

## 1. PubMed XML Format

### 1.1 DTD Structure

PubMed uses the \`pubmed_230101.dtd\` Document Type Definition (updated annually).

\`\`\`xml
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMed 230101//EN"
    "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_230101.dtd">
\`\`\`

### 1.2 Available Metadata Fields

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

PMC uses the Journal Article Tag Suite (JATS) standard for full-text content including:
- Complete article body with sections
- Tables, figures (with captions)
- Supplementary materials
- References with full citation data
- Author affiliations, Funding information

---

## 3. Core Data Structures (TypeScript)

### 3.1 Research Paper Interface

\`\`\`typescript
interface ResearchPaper {
  pmid: number;
  pmcid?: string;
  doi?: string;
  title: string;
  abstract: string;
  journal: Journal;
  publicationDate: Date;
  authors: Author[];
  meshTerms: MeSHTerm[];
  keywords: string[];
  mentionedGenes: string[];
  mentionedSNPs: string[];
  abstractEmbedding?: Float32Array;
  chunks?: PaperChunk[];
  citations: string[];
}
\`\`\`

---

## 4. RuVector Storage Schema

### 4.1 Collections

- **articles**: Main article metadata (dimension: 384, cosine distance)
- **paper_chunks**: Individual chunks for RAG (dimension: 384, cosine distance)

### 4.2 Graph Relationships

\`\`\`cypher
(:Article)-[:CITES]->(:Article)
(:Article)-[:MENTIONS_GENE]->(:Gene)
(:Article)-[:MENTIONS_SNP]->(:SNP)
(:PaperChunk)-[:CHUNK_OF]->(:Article)
\`\`\`

---

## 5. Size Estimates

### Database Size (5M Papers)

| Component | Size |
|-----------|------|
| Articles (metadata + abstract) | ~10 GB |
| Article embeddings (compressed) | ~800 MB |
| Chunks | ~10 GB |
| Chunk embeddings (compressed) | ~4 GB |
| Graph edges | ~5 GB |
| **Total** | **~30 GB** |

---

## Download

| Source | Method | URL/Command |
|--------|--------|-------------|
| PubMed XML | FTP | ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/ |
| PMC JATS | FTP | ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/ |

---

## License

| Source | License | Commercial Use |
|--------|---------|----------------|
| PubMed | Public domain | Yes |
| PMC OA | CC BY/CC0 | Yes (OA subset) |

---

## Glossary

| Term | Definition |
|------|------------|
| JATS XML | Journal Article Tag Suite - NISO standard for full-text |
| PMID | PubMed Identifier - unique integer for each citation |
| MeSH | Medical Subject Headings - controlled vocabulary |
| RAG | Retrieval-Augmented Generation |
| HNSW | Hierarchical Navigable Small World - vector index |

---

*Full content preserved from original source at databases/literature/data-structures.md*
