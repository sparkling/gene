/**
 * RuVector Graph Library for Gene Project
 *
 * Provides a TypeScript API for the THREE WORLDS knowledge graph:
 * - World 1: Genetics (Gene, SNP, Pathway, Condition)
 * - World 2: Traditional Medicine (Herb, Formula, Modality)
 * - World 3: Nutrition (Nutrient, Food, Biomarker)
 *
 * Features:
 * - Native @ruvector/graph-node bindings (10x faster than WASM)
 * - Cypher query language support
 * - Hypergraph support for multi-node relationships
 * - Vector similarity search with HNSW indexing
 * - ACID transactions
 * - Integration with ruvector vector embeddings
 *
 * @module ruvector-graph
 * @version 1.0.0
 */

// ============================================
// Type Definitions
// ============================================

/**
 * Distance metrics supported by RuVector
 */
export type DistanceMetric = 'Cosine' | 'Euclidean' | 'DotProduct';

/**
 * Graph database configuration
 */
export interface GraphConfig {
  /** Path for persistent storage (omit for in-memory) */
  storagePath?: string;
  /** Embedding dimensions (default: 384 for all-MiniLM-L6-v2) */
  dimensions?: number;
  /** Distance metric for similarity search */
  distanceMetric?: DistanceMetric;
}

/**
 * Node types in THREE WORLDS schema
 */
export type NodeLabel =
  // World 1: Genetics
  | 'Gene'
  | 'SNP'
  | 'Pathway'
  | 'Condition'
  | 'Publication'
  // World 2: Traditional Medicine
  | 'Herb'
  | 'Formula'
  | 'Modality'
  | 'TCMPattern'
  | 'Meridian'
  // World 3: Nutrition
  | 'Nutrient'
  | 'Food'
  | 'Biomarker'
  | 'Supplement'
  | 'FoodGroup';

/**
 * Edge types for relationships
 */
export type EdgeType =
  // Genetics internal
  | 'HAS_VARIANT'
  | 'INVOLVED_IN'
  | 'ASSOCIATED_WITH'
  | 'REGULATES'
  // Traditional medicine internal
  | 'INGREDIENT_OF'
  | 'BELONGS_TO'
  | 'ENTERS'
  | 'PATTERN_OF'
  // Nutrition internal
  | 'CONTAINS'
  | 'MEASURES'
  | 'IN_GROUP'
  // Cross-world relationships
  | 'AFFECTS'
  | 'TREATS'
  | 'REQUIRES_COFACTOR'
  | 'TARGETS'
  | 'MODULATES'
  | 'DEFICIENCY_CAUSES'
  | 'HAS_RESEARCH'
  | 'INTERACTS_WITH';

/**
 * Input for creating a node
 */
export interface NodeInput {
  /** Unique node identifier */
  id: string;
  /** 384-dimensional embedding vector */
  embedding: Float32Array | number[];
  /** Node labels (e.g., ['Gene'], ['SNP', 'Variant']) */
  labels: NodeLabel[];
  /** Node properties as key-value pairs */
  properties?: Record<string, string>;
}

/**
 * Input for creating an edge
 */
export interface EdgeInput {
  /** Source node ID */
  from: string;
  /** Target node ID */
  to: string;
  /** Edge type/description */
  description: EdgeType | string;
  /** Edge embedding for similarity search */
  embedding: Float32Array | number[];
  /** Confidence score (0-1) */
  confidence?: number;
  /** Additional metadata */
  metadata?: Record<string, unknown>;
}

/**
 * Input for creating a hyperedge (multi-node relationship)
 */
export interface HyperedgeInput {
  /** Node IDs participating in the hyperedge */
  nodes: string[];
  /** Description of the relationship */
  description: string;
  /** Hyperedge embedding */
  embedding: Float32Array | number[];
  /** Confidence score (0-1) */
  confidence?: number;
  /** Additional metadata */
  metadata?: Record<string, unknown>;
}

/**
 * Result from a Cypher query
 */
export interface QueryResult {
  nodes: Array<{
    id: string;
    labels: string[];
    properties: Record<string, string>;
  }>;
  edges: Array<{
    from: string;
    to: string;
    description: string;
    confidence: number;
  }>;
}

/**
 * Result from vector similarity search
 */
export interface SimilarityResult {
  id: string;
  score: number;
  description?: string;
  metadata?: Record<string, unknown>;
}

/**
 * Graph database statistics
 */
export interface GraphStats {
  totalNodes: number;
  totalEdges: number;
  totalHyperedges: number;
  avgDegree: number;
}

// ============================================
// RuVector Graph Database Class
// ============================================

/**
 * RuVector Graph Database wrapper for the Gene project
 *
 * @example
 * ```typescript
 * import { RuVectorGraph } from './lib/ruvector-graph';
 *
 * // Create persistent database
 * const graph = new RuVectorGraph({
 *   storagePath: './data/gene-knowledge.db',
 *   dimensions: 384,
 *   distanceMetric: 'Cosine'
 * });
 *
 * // Create a Gene node
 * await graph.createNode({
 *   id: 'gene_MTHFR',
 *   embedding: geneEmbedding,
 *   labels: ['Gene'],
 *   properties: { symbol: 'MTHFR', name: 'Methylenetetrahydrofolate reductase' }
 * });
 *
 * // Query with Cypher
 * const result = await graph.query('MATCH (g:Gene) RETURN g');
 * ```
 */
export class RuVectorGraph {
  private db: InstanceType<typeof GraphDatabase>;
  private config: Required<GraphConfig>;
  private initialized: boolean = false;

  constructor(config: GraphConfig = {}) {
    this.config = {
      storagePath: config.storagePath || '',
      dimensions: config.dimensions || 384,
      distanceMetric: config.distanceMetric || 'Cosine',
    };
  }

  /**
   * Initialize the graph database
   * Must be called before any operations
   */
  async initialize(): Promise<void> {
    if (this.initialized) {
      return;
    }

    const { GraphDatabase } = await import('@ruvector/graph-node');

    const options: Record<string, unknown> = {
      distanceMetric: this.config.distanceMetric,
      dimensions: this.config.dimensions,
    };

    if (this.config.storagePath) {
      options.storagePath = this.config.storagePath;
    }

    this.db = new GraphDatabase(options);
    this.initialized = true;
  }

  /**
   * Ensure database is initialized
   */
  private ensureInitialized(): void {
    if (!this.initialized) {
      throw new Error('Database not initialized. Call initialize() first.');
    }
  }

  // ==========================================
  // Node Operations
  // ==========================================

  /**
   * Create a node in the graph
   */
  async createNode(input: NodeInput): Promise<string> {
    this.ensureInitialized();

    const embedding = input.embedding instanceof Float32Array
      ? input.embedding
      : new Float32Array(input.embedding);

    if (embedding.length !== this.config.dimensions) {
      throw new Error(
        `Embedding must have ${this.config.dimensions} dimensions, got ${embedding.length}`
      );
    }

    return await this.db.createNode({
      id: input.id,
      embedding,
      labels: input.labels,
      properties: input.properties || {},
    });
  }

  /**
   * Batch create multiple nodes (131K+ ops/sec)
   */
  async batchCreateNodes(nodes: NodeInput[]): Promise<string[]> {
    this.ensureInitialized();

    const preparedNodes = nodes.map((node) => ({
      id: node.id,
      embedding: node.embedding instanceof Float32Array
        ? node.embedding
        : new Float32Array(node.embedding),
      labels: node.labels,
      properties: node.properties || {},
    }));

    const result = await this.db.batchInsert({
      nodes: preparedNodes,
      edges: [],
    });

    return result.nodeIds;
  }

  // ==========================================
  // Edge Operations
  // ==========================================

  /**
   * Create an edge between two nodes
   */
  async createEdge(input: EdgeInput): Promise<string> {
    this.ensureInitialized();

    const embedding = input.embedding instanceof Float32Array
      ? input.embedding
      : new Float32Array(input.embedding);

    return await this.db.createEdge({
      from: input.from,
      to: input.to,
      description: input.description,
      embedding,
      confidence: input.confidence ?? 1.0,
      metadata: input.metadata || {},
    });
  }

  /**
   * Batch create multiple edges
   */
  async batchCreateEdges(edges: EdgeInput[]): Promise<string[]> {
    this.ensureInitialized();

    const preparedEdges = edges.map((edge) => ({
      from: edge.from,
      to: edge.to,
      description: edge.description,
      embedding: edge.embedding instanceof Float32Array
        ? edge.embedding
        : new Float32Array(edge.embedding),
      confidence: edge.confidence ?? 1.0,
      metadata: edge.metadata || {},
    }));

    const result = await this.db.batchInsert({
      nodes: [],
      edges: preparedEdges,
    });

    return result.edgeIds;
  }

  // ==========================================
  // Hyperedge Operations (Multi-node relationships)
  // ==========================================

  /**
   * Create a hyperedge connecting multiple nodes
   *
   * Use cases:
   * - Formula compositions (multiple herbs)
   * - Gene-Nutrient-Condition pathways
   * - Complex pathway interactions
   */
  async createHyperedge(input: HyperedgeInput): Promise<string> {
    this.ensureInitialized();

    if (input.nodes.length < 2) {
      throw new Error('Hyperedge requires at least 2 nodes');
    }

    const embedding = input.embedding instanceof Float32Array
      ? input.embedding
      : new Float32Array(input.embedding);

    return await this.db.createHyperedge({
      nodes: input.nodes,
      description: input.description,
      embedding,
      confidence: input.confidence ?? 1.0,
      metadata: input.metadata || {},
    });
  }

  /**
   * Search hyperedges by vector similarity
   */
  async searchHyperedges(
    embedding: Float32Array | number[],
    k: number = 10
  ): Promise<SimilarityResult[]> {
    this.ensureInitialized();

    const embeddingArray = embedding instanceof Float32Array
      ? embedding
      : new Float32Array(embedding);

    return await this.db.searchHyperedges({
      embedding: embeddingArray,
      k,
    });
  }

  // ==========================================
  // Query Operations
  // ==========================================

  /**
   * Execute a Cypher query
   *
   * @example
   * ```typescript
   * // Find all genes
   * const genes = await graph.query('MATCH (g:Gene) RETURN g');
   *
   * // Find SNPs affecting nutrients
   * const snpNutrients = await graph.query(`
   *   MATCH (s:SNP)-[a:AFFECTS]->(n:Nutrient)
   *   RETURN s, a, n
   * `);
   *
   * // Find herbs treating conditions via pathways
   * const herbConditions = await graph.query(`
   *   MATCH (h:Herb)-[:TARGETS]->(p:Pathway)<-[:INVOLVED_IN]-(g:Gene)-[:HAS_VARIANT]->(s:SNP)-[:ASSOCIATED_WITH]->(c:Condition)
   *   RETURN h, p, g, s, c
   * `);
   * ```
   */
  async query(cypher: string): Promise<QueryResult> {
    this.ensureInitialized();
    return await this.db.query(cypher);
  }

  /**
   * Execute a Cypher query synchronously
   */
  querySync(cypher: string): QueryResult {
    this.ensureInitialized();
    return this.db.querySync(cypher);
  }

  // ==========================================
  // Graph Traversal
  // ==========================================

  /**
   * Get k-hop neighbors from a starting node
   *
   * @example
   * ```typescript
   * // Get direct neighbors of MTHFR gene
   * const neighbors1 = await graph.getNeighbors('gene_MTHFR', 1);
   *
   * // Get 2-hop neighborhood
   * const neighbors2 = await graph.getNeighbors('gene_MTHFR', 2);
   * ```
   */
  async getNeighbors(nodeId: string, k: number = 1): Promise<string[]> {
    this.ensureInitialized();
    return await this.db.kHopNeighbors(nodeId, k);
  }

  // ==========================================
  // Transaction Operations
  // ==========================================

  /**
   * Begin a transaction
   */
  async beginTransaction(): Promise<string> {
    this.ensureInitialized();
    return await this.db.begin();
  }

  /**
   * Commit a transaction
   */
  async commitTransaction(txId: string): Promise<void> {
    this.ensureInitialized();
    await this.db.commit(txId);
  }

  /**
   * Rollback a transaction
   */
  async rollbackTransaction(txId: string): Promise<void> {
    this.ensureInitialized();
    await this.db.rollback(txId);
  }

  /**
   * Execute operations within a transaction
   */
  async withTransaction<T>(
    fn: (txId: string) => Promise<T>
  ): Promise<T> {
    const txId = await this.beginTransaction();
    try {
      const result = await fn(txId);
      await this.commitTransaction(txId);
      return result;
    } catch (error) {
      await this.rollbackTransaction(txId);
      throw error;
    }
  }

  // ==========================================
  // Statistics & Utilities
  // ==========================================

  /**
   * Get graph statistics
   */
  async getStats(): Promise<GraphStats> {
    this.ensureInitialized();
    return await this.db.stats();
  }

  /**
   * Check if database is persistent
   */
  isPersistent(): boolean {
    this.ensureInitialized();
    return this.db.isPersistent();
  }

  /**
   * Get storage path
   */
  getStoragePath(): string | null {
    this.ensureInitialized();
    return this.db.getStoragePath();
  }

  /**
   * Get configuration
   */
  getConfig(): Required<GraphConfig> {
    return { ...this.config };
  }
}

// ============================================
// Embedding Utilities
// ============================================

/**
 * Generate a zero embedding (placeholder)
 */
export function zeroEmbedding(dimensions: number = 384): Float32Array {
  return new Float32Array(dimensions);
}

/**
 * Generate a normalized random embedding (for testing)
 */
export function randomEmbedding(dimensions: number = 384): Float32Array {
  const embedding = new Float32Array(dimensions);
  let norm = 0;

  for (let i = 0; i < dimensions; i++) {
    embedding[i] = (Math.random() - 0.5) * 2;
    norm += embedding[i] * embedding[i];
  }

  norm = Math.sqrt(norm);
  for (let i = 0; i < dimensions; i++) {
    embedding[i] /= norm;
  }

  return embedding;
}

/**
 * Normalize an embedding to unit length
 */
export function normalizeEmbedding(embedding: Float32Array): Float32Array {
  let norm = 0;
  for (let i = 0; i < embedding.length; i++) {
    norm += embedding[i] * embedding[i];
  }
  norm = Math.sqrt(norm);

  const normalized = new Float32Array(embedding.length);
  for (let i = 0; i < embedding.length; i++) {
    normalized[i] = embedding[i] / norm;
  }

  return normalized;
}

/**
 * Calculate cosine similarity between two embeddings
 */
export function cosineSimilarity(a: Float32Array, b: Float32Array): number {
  if (a.length !== b.length) {
    throw new Error('Embeddings must have same dimensions');
  }

  let dot = 0;
  let normA = 0;
  let normB = 0;

  for (let i = 0; i < a.length; i++) {
    dot += a[i] * b[i];
    normA += a[i] * a[i];
    normB += b[i] * b[i];
  }

  return dot / (Math.sqrt(normA) * Math.sqrt(normB));
}

// ============================================
// GraphDatabase Type Declaration (for dynamic import)
// ============================================

declare const GraphDatabase: new (options: Record<string, unknown>) => {
  createNode: (input: unknown) => Promise<string>;
  createEdge: (input: unknown) => Promise<string>;
  createHyperedge: (input: unknown) => Promise<string>;
  batchInsert: (input: unknown) => Promise<{ nodeIds: string[]; edgeIds: string[] }>;
  searchHyperedges: (input: unknown) => Promise<SimilarityResult[]>;
  query: (cypher: string) => Promise<QueryResult>;
  querySync: (cypher: string) => QueryResult;
  kHopNeighbors: (nodeId: string, k: number) => Promise<string[]>;
  begin: () => Promise<string>;
  commit: (txId: string) => Promise<void>;
  rollback: (txId: string) => Promise<void>;
  stats: () => Promise<GraphStats>;
  isPersistent: () => boolean;
  getStoragePath: () => string | null;
};

// ============================================
// Default Export
// ============================================

export default RuVectorGraph;
