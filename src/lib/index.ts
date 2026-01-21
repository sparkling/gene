/**
 * RuVector Graph Library for Gene Project
 *
 * This module provides a complete TypeScript API for building and querying
 * the THREE WORLDS knowledge graph:
 *
 * - World 1: Genetics (Gene, SNP, Pathway, Condition)
 * - World 2: Traditional Medicine (Herb, Formula, Modality)
 * - World 3: Nutrition (Nutrient, Food, Biomarker)
 *
 * @module lib
 * @version 1.0.0
 *
 * @example
 * ```typescript
 * import {
 *   RuVectorGraph,
 *   RuVectorEmbeddings,
 *   randomEmbedding
 * } from './lib';
 *
 * // Initialize graph and embeddings
 * const graph = new RuVectorGraph({ storagePath: './data/knowledge.db' });
 * const embeddings = new RuVectorEmbeddings();
 *
 * await graph.initialize();
 * await embeddings.initialize();
 *
 * // Create a gene node
 * const geneEmbed = await embeddings.embedGene({
 *   symbol: 'MTHFR',
 *   name: 'Methylenetetrahydrofolate reductase',
 *   description: 'Key enzyme in folate metabolism'
 * });
 *
 * await graph.createNode({
 *   id: 'gene_MTHFR',
 *   embedding: geneEmbed,
 *   labels: ['Gene'],
 *   properties: { symbol: 'MTHFR' }
 * });
 *
 * // Query with Cypher
 * const result = await graph.query('MATCH (g:Gene) RETURN g');
 * ```
 */

// ============================================
// Graph Database
// ============================================

export {
  RuVectorGraph,
  // Types
  type DistanceMetric,
  type GraphConfig,
  type NodeLabel,
  type EdgeType,
  type NodeInput,
  type EdgeInput,
  type HyperedgeInput,
  type QueryResult,
  type SimilarityResult,
  type GraphStats,
  // Utilities
  zeroEmbedding,
  randomEmbedding,
  normalizeEmbedding,
  cosineSimilarity,
} from './ruvector-graph';

// ============================================
// Embeddings
// ============================================

export {
  RuVectorEmbeddings,
  createEmbeddings,
  // Types
  type EmbeddingModel,
  type EmbeddingConfig,
  type HNSWConfig,
  type VectorSearchResult,
  type HyperbolicConfig,
} from './ruvector-embeddings';

// ============================================
// Demo
// ============================================

export { runDemo } from './ruvector-graph-demo';

// ============================================
// Default Export
// ============================================

import RuVectorGraph from './ruvector-graph';
export default RuVectorGraph;
