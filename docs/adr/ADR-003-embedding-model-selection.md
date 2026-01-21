# ADR-003: Embedding Model Selection

## Status
Accepted

## Context

The Gene platform needs vector embeddings for:
1. **Gene similarity search** - Find genes with similar functions
2. **SNP similarity search** - Find SNPs with similar clinical significance
3. **AI agent memory** - Claude Flow patterns and trajectories
4. **RAG retrieval** - Document and knowledge retrieval

Key challenge: Gene symbols (MTHFR) and SNP IDs (rs1801133) are bare identifiers with no semantic meaning when embedded alone.

## Decision

### Model Selection: Start Simple

Use **all-mpnet-base-v2** (768 dimensions) as the primary embedding model.

| Phase | Model | Dimensions | Use For |
|-------|-------|------------|---------|
| **Phase 1 (Now)** | all-mpnet-base-v2 | 768 | Everything - solid for biomedical text |
| **Phase 2 (If needed)** | Add PubMedBERT | 768 | Better gene/SNP description similarity |
| **Phase 3 (If needed)** | Add Gene2Vec | 200 | Gene-gene functional relationships |
| **Phase 4 (If sequences)** | Add ESM-2 | 320-480 | Protein sequence similarity |

### Key Principle: Always Embed with Context

Gene symbols and SNP IDs alone don't embed well. Always enrich with context:

```typescript
// ❌ POOR - bare identifier has no semantic meaning
await embed("MTHFR");

// ✅ GOOD - context-enriched embedding
await embed("MTHFR: Methylenetetrahydrofolate reductase, key enzyme in methylation cycle");

// ✅ GOOD - SNP with clinical significance
await embed("rs1801133 (MTHFR C677T): Reduces enzyme activity, impacts folate requirements");
```

### Embedding Implementation

```typescript
// Enrich gene data before embedding
function enrichGeneText(gene: Gene): string {
  return `${gene.symbol}: ${gene.name}. ${gene.description || ''}. ` +
         `Pathway: ${gene.pathway || 'unknown'}`;
}

// Enrich SNP data before embedding
function enrichSNPText(snp: SNP): string {
  return `${snp.rsid} (${snp.commonName || snp.rsid}): ` +
         `${snp.clinicalSignificance || 'unknown significance'}. ` +
         `Gene: ${snp.geneSymbol || 'unknown'}`;
}
```

## Consequences

### Positive
- Simple architecture (one model, one embedding space)
- all-mpnet-base-v2 is FREE and open source (Sentence Transformers)
- Research shows general models often match specialized biomedical models
- Already proven in RuVector stack
- Can add specialized models later in separate namespaces

### Negative
- May not capture gene-gene functional relationships as well as Gene2Vec
- May not capture protein structure similarity as well as ESM-2
- Requires context enrichment for good results

### Neutral
- 768 dimensions is standard size
- Need to maintain context enrichment functions

## Options Considered

### Option 1: all-mpnet-base-v2 Only (Chosen for Phase 1)
- **Pros**: Simple, proven, FREE, good general performance
- **Cons**: Not specialized for biomedical domain

### Option 2: PubMedBERT Only
- **Pros**: Trained on biomedical literature
- **Cons**: May not generalize as well, more complex setup

### Option 3: Multiple Specialized Models
- **Pros**: Best-in-class for each data type
- **Cons**: Complex architecture, multiple embedding spaces, harder to maintain

### Option 4: OpenAI text-embedding-3-small
- **Pros**: High quality, easy API
- **Cons**: API costs, vendor dependency, not self-hosted

## Future Expansion (If Needed)

If similarity results aren't good enough, add specialized models:

| Model | When to Add | What It Improves |
|-------|-------------|------------------|
| **PubMedBERT** | Gene descriptions not clustering well | Biomedical text understanding |
| **Gene2Vec** | Need gene-gene functional similarity | Co-expression relationships |
| **SNP2Vec** | SNP risk prediction accuracy low | Variant-specific embeddings |
| **ESM-2** | Have protein sequences to compare | Amino acid sequence similarity |

## Related Decisions
- [ADR-001](./ADR-001-three-database-architecture.md): Three-Database Architecture
- [ADR-005](./ADR-005-claude-flow-memory.md): Claude Flow Memory Configuration

## References
- [Sentence Transformers](https://www.sbert.net/) - all-mpnet-base-v2
- [PubMedBERT](https://huggingface.co/microsoft/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext)
- [Gene2Vec](https://www.nature.com/articles/s41598-018-23495-9)
- [ESM-2](https://github.com/facebookresearch/esm) - Protein language model
