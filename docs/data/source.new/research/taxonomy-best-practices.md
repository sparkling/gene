---
id: taxonomy-best-practices
title: Taxonomy Best Practices for Scientific Data Source Classification
category: research
tier: 1
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [taxonomy, classification, ontology, polyhierarchy, faceted-classification, best-practices]
---

# Taxonomy Best Practices for Scientific Data Source Classification

> **Last Updated**: January 2026
> **Purpose**: Research-backed guidelines for building effective taxonomies to categorize scientific and biomedical data sources.
> **Audience**: Queen Taxonomy Architect and taxonomy design team

---

## Executive Summary

This document synthesizes best practices from information science, biomedical ontologies, and taxonomy research to guide the design of a polyhierarchical taxonomy for scientific data sources. Key findings:

| Principle | Recommendation |
|-----------|---------------|
| **Structure** | Polyhierarchical with 3-7 facets |
| **Depth** | 2-4 levels per facet (max) |
| **Top-level categories** | 10-15 maximum |
| **Approach** | Hybrid (top-down + bottom-up) |
| **MECE adaptation** | Collectively Exhaustive is higher priority than Mutually Exclusive |

---

## 1. Core Taxonomy Principles

### 1.1 Polyhierarchical vs. Monohierarchical Structures

**Definition**: Polyhierarchy occurs when a concept has two or more broader concepts within a taxonomy hierarchy. This allows different users to navigate from different starting points to reach the same concept.

**When to Use Polyhierarchy** (from [Hedden Information Management](https://www.hedden-information.com/polyhierarchy-in-taxonomies/)):
- When "the concept's relationship is correctly and inherently hierarchical in both of its cases"
- Example: "Educational software" legitimately belongs under both "Software" AND "Educational products"

**Best Practice Guidelines**:
- Polyhierarchies typically involve only 2 broader concepts
- Multiple instances of 3+ parents suggest design problems
- In faceted taxonomies, limit polyhierarchies to 2-3 maximum across the entire system
- Use the "All-and-Some" test: if "All X are Y" AND "Some Y are X" holds true for both parent relationships, polyhierarchy is valid

**For Data Sources**: A database like "ClinVar" could legitimately appear under:
- "Variant/Mutation Databases"
- "Clinical Databases"
- "NCBI Resources"

This polyhierarchy serves users who approach from different mental models.

### 1.2 Adapting MECE for Overlapping Categories

**Traditional MECE** (Mutually Exclusive, Collectively Exhaustive) requires:
- Each item appears in exactly one category
- All items have a category

**Adapted MECE for Scientific Data Sources**:

According to research on the [MECE principle](https://en.wikipedia.org/wiki/MECE_principle), "if you had to choose one over the other, being Collectively Exhaustive should be a higher priority than being Mutually Exclusive."

**Practical Application**:

| Principle | Strict MECE | Adapted for Data Sources |
|-----------|-------------|--------------------------|
| Mutually Exclusive | Each source in ONE category | Sources can appear in multiple FACETS |
| Collectively Exhaustive | All sources covered | Maintained - all sources must have at least one category |
| Overlap handling | Forbidden | Managed through faceted structure |

**Resolution Strategy**:
1. **Within a single facet**: Maintain mutual exclusivity (a source is ONE data type)
2. **Across facets**: Allow overlap (a source can be both "Genomic" AND "Clinical")
3. **Use facets for different classification dimensions** rather than forcing single-category placement

### 1.3 Granularity Decisions

**Key Insight** from [Hedden Information Management](https://www.hedden-information.com/taxonomy-hierarchy-levels/): "Taxonomies naturally have hierarchies, but do not naturally have levels, which are an artificial layer."

**When to Stop Subdividing**:

1. **The One-Child Rule**: A term should have no narrower terms OR at least two narrower terms, but never just one
2. **User Value Test**: Does the subdivision help users find what they need faster?
3. **Content Threshold**: A category should have meaningful content assigned (avoid "singularities")
4. **Navigation Efficiency**: "If users need more than three clicks to find common information, restructure top-level categories rather than adding hierarchy depth"

**Recommended Depth**:

| Context | Recommended Depth |
|---------|-------------------|
| Top-level facets | 10-15 maximum |
| Each facet hierarchy | 2-4 levels deep |
| Total clicks to content | 3 maximum |

**Granularity Types** (from [Keet's taxonomy of granularity](http://www.meteck.org/files/GrCGranTypes_CMK.pdf)):
- **Spatial**: Geographic or structural decomposition
- **Conceptual**: Abstraction levels (general to specific)
- **Temporal**: Time-based divisions
- **Functional**: Purpose-based groupings

### 1.4 Faceted Classification Systems

**Definition**: A faceted taxonomy consists of multiple independent hierarchies (facets) that describe different aspects of the same resource.

**Core Principles** (from [ScienceDirect](https://www.sciencedirect.com/topics/computer-science/faceted-classification) and [Hedden](https://www.hedden-information.com/faceted-classification-and-faceted-taxonomies/)):

1. **Mutual exclusivity within facets**: Each concept belongs to only ONE facet
2. **Independence across facets**: Facets can be combined freely
3. **Optimal facet count**: 3-7 facets is practical and manageable
4. **Question-based design**: Each facet answers a different question (What? Who? Where? When?)

**Proposed Facets for Data Sources**:

| Facet | Question Answered | Example Values |
|-------|-------------------|----------------|
| **Data Type** | What kind of data? | Sequence, Variant, Expression, Structure |
| **Domain** | What biological domain? | Genomics, Proteomics, Metabolomics, Clinical |
| **Application** | What is it used for? | Research, Clinical, Drug Discovery |
| **Organism** | What species? | Human, Model Organisms, Multi-species |
| **Curation Level** | How is data quality assured? | Expert-curated, Community-curated, Automated |
| **Access Model** | How to get data? | Open, Registered, Controlled |

---

## 2. Taxonomy Design Patterns

### 2.1 Top-Down vs. Bottom-Up Approaches

**Top-Down Approach** (from [ANSI/NISO Z39.19](https://www.niso.org/standards-committees/vocabulary-standards)):
- Start with broadest terms
- Define narrower terms to reach desired specificity
- Build hierarchical relationships as work proceeds

**Strengths**:
- Ensures conceptual coherence
- Aligns with organizational strategy
- Good for well-understood domains

**Weaknesses**:
- May miss important concepts that emerge from actual content
- Risk of imposing structure that doesn't match reality

**Bottom-Up Approach**:
- Derive terms from actual content corpus
- Group related terms into broader categories
- Build hierarchy from specific to general

**Strengths**:
- Grounded in real content
- Captures user language
- Discovers unexpected categories

**Weaknesses**:
- May lack conceptual coherence
- Difficult to manage at scale

### 2.2 Hybrid Approach (Recommended)

According to [Enterprise Knowledge](https://enterprise-knowledge.com/the-art-of-taxonomy-design/), "the hybrid taxonomy design methodology leverages both top-down engagement and bottom-up analyses."

**Process**:

```
Phase 1: Top-Down Foundation
├── Review existing scientific ontologies (MeSH, GO, OBO)
├── Identify major data type categories from literature
├── Conduct stakeholder interviews
└── Define initial facet structure

Phase 2: Bottom-Up Validation
├── Analyze actual data sources in inventory
├── Extract terminology from source documentation
├── Review user search logs and queries
└── Identify gaps and emergent categories

Phase 3: Integration
├── Reconcile top-down structure with bottom-up findings
├── Validate with domain experts
├── Test with representative users
└── Iterate based on feedback
```

### 2.3 Determining Category Boundaries

**Questions to Ask**:
1. Do users think of these as the same or different?
2. Do the items share the same attributes?
3. Would combining them help or hurt findability?
4. Is the distinction meaningful for the use case?

**Boundary Tests** (adapted from [NN/Group](https://www.nngroup.com/articles/taxonomy-101/)):

| Test | Pass Criteria |
|------|---------------|
| **Card Sort** | 80%+ agreement on category placement |
| **Find-It Test** | Users locate items in <3 attempts |
| **Name Test** | Category name immediately understood |
| **Sibling Test** | Categories at same level are peers |

### 2.4 Naming Conventions

**Best Practices** (from [NN/Group](https://www.nngroup.com/articles/taxonomy-101/)):

1. **Use user language**: Match how users actually search and talk
2. **Be specific**: "Variant Annotation Databases" > "Variant Data" > "Variants"
3. **Be consistent**: If using plurals, use everywhere; same for capitalization
4. **Avoid jargon**: Or provide clear definitions
5. **Include synonyms**: As non-preferred terms pointing to preferred terms

**For Scientific Data Sources**:
- Use established terminology from MeSH, Gene Ontology where applicable
- Include common abbreviations as synonyms (e.g., "GO" -> "Gene Ontology")
- Prefer descriptive names that indicate content (e.g., "Protein-Protein Interaction Databases")

---

## 3. Common Pitfalls to Avoid

### 3.1 Over-Categorization

**Problem**: Creating too many categories or too deep hierarchies.

**Signs of Over-Categorization** (from [Earley Information Science](https://www.earley.com/insights/ten-common-mistakes-when-developing-taxonomy)):
- Categories with only 1-2 items
- Users can't find items
- Taggers disagree on placement
- Deep hierarchies require many clicks

**Solution**: Every category must have a documented use case. If a distinction doesn't help users, collapse the categories.

### 3.2 Category Overlap Issues

**Problem**: Ambiguous items that fit multiple categories, causing inconsistent tagging.

**Signs** (from [Matrix Flows](https://www.matrixflows.com/blog/10-best-practices-for-creating-taxonomy-for-your-company-knowledge-base)):
- Cross-listing appears frequently (>20% of items)
- Same item tagged differently by different people
- Users search multiple paths for same content

**Solutions**:
1. Use faceted structure to handle legitimate multi-dimensionality
2. Limit polyhierarchy to 2 parent concepts maximum
3. Write clear scope notes for boundary cases
4. Create decision trees for ambiguous cases

### 3.3 Orphan Categories

**Problem**: Categories with no content or very few items.

**Causes**:
- Speculative categories created "in case we need them"
- Overly granular subdivision
- Content that no longer exists

**Solution**: Apply content threshold - if a category has <3 items, merge it with siblings or parent.

### 3.4 Imbalanced Trees

**Problem**: Some branches very deep/full, others shallow/sparse.

**Expected Distribution** (from [Taxonomy Strategies](https://taxonomystrategies.com/taxonomy-assessment-quantitative-analytics-practices-heuristics/)):
Content naturally follows a **Zipf distribution**:
- Most frequent category: 2x content of second
- Second: 3x content of third
- And so on

**When to Worry**:
- A category has >50% of all content -> Split it
- A category has <1% of content -> Merge or delete it
- Sibling categories have wildly different depths -> Rebalance

### 3.5 Additional Pitfalls (from [Earley](https://www.earley.com/insights/ten-common-mistakes-when-developing-taxonomy))

| Pitfall | Description | Mitigation |
|---------|-------------|------------|
| **Confusing taxonomy with navigation** | Navigation is one use; taxonomy enables multiple views | Design for multiple use cases |
| **Adopting generic taxonomies** | Pre-built taxonomies miss organizational context | Customize to your content and users |
| **Neglecting maintenance** | Taxonomies require ongoing governance | Establish review cadence and ownership |
| **Poor implementation** | Taxonomy "bolted on" rather than integrated | Build into core architecture |
| **Insufficient tagging compliance** | Users don't apply taxonomy consistently | Train users, provide guidance |

---

## 4. Biomedical/Scientific Taxonomy Examples

### 4.1 Gene Ontology (GO)

**Structure** (from [Gene Ontology Resource](https://geneontology.org/)):
- Three main facets (aspects):
  - Molecular Function
  - Biological Process
  - Cellular Component
- ~45,000 terms
- Polyhierarchical: terms can have multiple parents

**What Makes It Successful**:
- Clear scope definitions for each term
- Formal relations between terms (is_a, part_of, regulates)
- Community curation and expert review
- Standardized identifiers (GO:0000000)
- Multiple formats (OWL, OBO, JSON)

### 4.2 OBO Foundry Principles

**Key Principles** (from [OBO Foundry](http://obofoundry.org/) and [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2814061/)):

1. **Openness**: Free to use without license constraints
2. **Orthogonality**: Each ontology has distinct, non-overlapping domain
3. **Unique identifiers**: Cross-ontology compatibility
4. **Formal definitions**: Genus-species definitions
5. **Versioning**: Track changes over time
6. **Bounded scope**: Don't mix unrelated domains

**Cross-Ontology Coordination**:
- Use Relation Ontology (RO) for standard relationships
- Build compound terms from existing ontology terms
- Maintain backward compatibility with identifiers

### 4.3 MeSH (Medical Subject Headings)

**Structure**:
- 16 top-level categories (Anatomy, Organisms, Diseases, etc.)
- ~30,000 descriptors
- Hierarchical tree structure with polyhierarchy
- Entry terms (synonyms) pointing to preferred terms

**Useful Patterns for Data Sources**:
- Qualifier system for adding context
- Scope notes explaining term usage
- Related terms for discovery

### 4.4 NCBI Database Classification

**Current Categorization** (from [NCBI All Resources](https://www.ncbi.nlm.nih.gov/guide/all/)):

| Category | Examples |
|----------|----------|
| DNA & RNA | GenBank, RefSeq, SRA |
| Proteins | UniProt, PDB, RefSeq Protein |
| Genes | Gene, OMIM |
| Genomes | Genome, Assembly |
| Chemicals | PubChem |
| Taxonomy | Taxonomy Browser |
| Literature | PubMed, PMC |
| Clinical | ClinVar, GTR |

**Lessons**:
- Categories based on data type (molecule-centric)
- Cross-links between related resources
- Hierarchy within each category

### 4.5 NIH Repository Classification

**Categories** (from [Scientific Data](https://www.nature.com/articles/s41597-024-03449-z)):
- Domain-specific repositories (categorized by subject)
- Open vs. Registered vs. Controlled access
- Generalist vs. Specialist repositories

---

## 5. Evaluation Criteria

### 5.1 Taxonomy Quality Attributes

**Seven Core Attributes** (from [Unterkalmsteiner et al., Expert Systems 2023](https://onlinelibrary.wiley.com/doi/full/10.1111/exsy.13098)):

| Attribute | Definition | Measurement |
|-----------|------------|-------------|
| **Comprehensiveness** | Capacity to categorize all objects in domain | % of items successfully categorized |
| **Robustness** | Ability to distinguish between items | Orthogonality and cohesiveness of categories |
| **Conciseness** | Classification using minimal concepts | Categories actually used vs. total categories |
| **Extensibility** | Ability to accommodate changes | Ease of adding/removing categories |
| **Explanatory** | Allows inferring characteristics from placement | User comprehension tests |
| **Mutual Exclusiveness** | Objects don't exist under multiple categories | % of items in single category (within facet) |
| **Reliability** | Different people classify same way | Inter-rater agreement (Kappa, Krippendorff's alpha) |

**Priority for Classification Tasks**:
Robustness > Comprehensiveness > Mutual Exclusiveness > Conciseness > Explanatory > Extensibility

### 5.2 Structural Evaluation Metrics

**Generalized Metrics** (from [EASE 2022](https://fuchss.org/assets/pdf/2022/ease-22-taxonomies.pdf)):

| Metric | What It Measures |
|--------|------------------|
| **Laconicity** | Avoids unnecessary categories |
| **Lucidity** | Categories are clearly distinguishable |
| **Completeness** | Covers the full domain |
| **Soundness** | Categories are logically valid |

### 5.3 Balance Testing

**Content Distribution Analysis** (from [Taxonomy Strategies](https://taxonomystrategies.com/taxonomy-assessment-quantitative-analytics-practices-heuristics/)):

1. Track items per category
2. Compare to Zipf distribution expectation
3. Categories exceeding expected frequency -> candidates to split
4. Categories below expected frequency -> candidates to merge

**Imbalance Indicators**:
- Any category >30% of total items
- Any category <1% of total items (unless specialized)
- Depth variance >2 levels between sibling branches

### 5.4 User-Centric Evaluation

**Testing Methods**:

| Method | What It Tests | Success Criteria |
|--------|---------------|------------------|
| **Card Sorting** | Category boundaries | 80%+ agreement |
| **Tree Testing** | Navigability | 70%+ success rate |
| **First-Click Testing** | Initial findability | >60% first click correct |
| **Search Log Analysis** | Category value | Users find via taxonomy, not just search |

**Qualitative Assessment** (from [ISO 25964](https://www.iso.org/standard/53657.html)):

1. **Structure**: Grouping principles, mutual exclusivity, consistency, depth/breadth
2. **Labeling**: Formatting consistency, editorial guidelines applied
3. **Attributes**: Clear definitions, scope notes, relationships
4. **Implementation**: Integration quality, usability for tagging

---

## 6. Recommendations for Data Source Taxonomy

### 6.1 Proposed Facet Structure

Based on the research above, recommend a **6-facet polyhierarchical structure**:

```
Facet 1: DATA TYPE (What kind of data?)
├── Sequence
│   ├── Nucleotide
│   ├── Protein
│   └── Structural
├── Variant/Mutation
├── Expression
├── Interaction
├── Pathway
├── Clinical
├── Literature
└── Reference/Annotation

Facet 2: DOMAIN (What biological area?)
├── Genomics
├── Transcriptomics
├── Proteomics
├── Metabolomics
├── Pharmacology
├── Phenotype/Disease
└── Taxonomy/Evolution

Facet 3: APPLICATION (Primary use case?)
├── Basic Research
├── Clinical/Diagnostic
├── Drug Discovery
├── Biomarker Discovery
├── Functional Annotation
└── Knowledge Integration

Facet 4: ORGANISM SCOPE
├── Human-specific
├── Model Organisms
├── Pathogen
├── Multi-species
└── Pan-taxonomic

Facet 5: CURATION MODEL
├── Expert-curated
├── Community-curated
├── Automated/Computed
└── Hybrid

Facet 6: ACCESS MODEL
├── Open (no restrictions)
├── Registered (account required)
├── Controlled (application required)
└── Commercial
```

### 6.2 Design Principles to Adopt

1. **Use hybrid approach**: Start with top-down structure informed by OBO/MeSH, validate with actual data sources
2. **Limit polyhierarchy**: Max 2 parents within any single facet
3. **Allow cross-facet tagging**: Each source gets one tag per facet (facets are independent)
4. **Target 2-3 levels depth**: Within each facet
5. **Include scope notes**: For every category, define what belongs and what doesn't
6. **Plan for governance**: Establish review process and ownership

### 6.3 Validation Plan

1. **Expert Review**: Subject matter experts validate category definitions
2. **Card Sort**: Test category boundaries with representative users
3. **Coverage Check**: Ensure all existing data sources can be categorized
4. **Balance Check**: Verify no category is over/under-populated
5. **Inter-rater Reliability**: Multiple taggers classify same sources, measure agreement

---

## 7. Sources

### Primary Sources

- [Hedden Information Management - Polyhierarchy](https://www.hedden-information.com/polyhierarchy-in-taxonomies/)
- [Hedden Information Management - Faceted Classification](https://www.hedden-information.com/faceted-classification-and-faceted-taxonomies/)
- [Hedden Information Management - Hierarchy Levels](https://www.hedden-information.com/taxonomy-hierarchy-levels/)
- [NN/Group - Taxonomy 101](https://www.nngroup.com/articles/taxonomy-101/)
- [Taxonomy Strategies - Assessment](https://taxonomystrategies.com/taxonomy-assessment-quantitative-analytics-practices-heuristics/)
- [Enterprise Knowledge - Taxonomy Design](https://enterprise-knowledge.com/the-art-of-taxonomy-design/)
- [Earley Information Science - Ten Common Mistakes](https://www.earley.com/insights/ten-common-mistakes-when-developing-taxonomy)

### Biomedical Ontology Sources

- [Gene Ontology Resource](https://geneontology.org/)
- [OBO Foundry](http://obofoundry.org/)
- [OBO Foundry Principles - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2814061/)
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
- [NCBI All Resources](https://www.ncbi.nlm.nih.gov/guide/all/)

### Research Papers

- [Unterkalmsteiner et al. - Taxonomy Quality Attributes (2023)](https://onlinelibrary.wiley.com/doi/full/10.1111/exsy.13098)
- [Biomedical Data Repositories - Nature Scientific Data](https://www.nature.com/articles/s41597-024-03449-z)
- [Biomedical Data Categorization - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10149933/)
- [Health Information Systems Taxonomy - JMIR](https://www.jmir.org/2024/1/e47682)

### Standards

- [ISO 25964 - Thesauri and interoperability](https://www.iso.org/standard/53657.html)
- [ANSI/NISO Z39.19 - Controlled Vocabularies](https://www.niso.org/standards-committees/vocabulary-standards)

---

## Appendix: Quick Reference Checklist

### Before Building

- [ ] Define primary use cases
- [ ] Identify user groups and their mental models
- [ ] Review existing ontologies (MeSH, GO, OBO Foundry)
- [ ] Inventory content to be categorized
- [ ] Decide facet structure (3-7 facets recommended)

### During Building

- [ ] Use hybrid top-down/bottom-up approach
- [ ] Limit depth to 2-4 levels per facet
- [ ] Apply "one-child rule" (never single child categories)
- [ ] Write scope notes for each category
- [ ] Limit polyhierarchy to 2 parents max
- [ ] Create preferred term + synonyms
- [ ] Document relationships (is-a, part-of, related-to)

### After Building

- [ ] Test with card sort (target 80% agreement)
- [ ] Check content distribution (Zipf pattern)
- [ ] Measure inter-rater reliability (Kappa >0.6)
- [ ] Validate comprehensiveness (all items categorized)
- [ ] Plan governance and maintenance cadence
- [ ] Train taggers and establish guidelines
