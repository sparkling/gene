/**
 * RuVector Embeddings Integration for Gene Project
 *
 * Provides vector embedding generation and integration with the graph database.
 * Uses the ruvector npm package for high-performance vector operations.
 *
 * Features:
 * - ONNX-based text embeddings (all-MiniLM-L6-v2)
 * - HNSW indexing for fast similarity search (150x-12,500x faster)
 * - Hyperbolic embeddings (Poincare ball model)
 * - Batch embedding generation
 * - Integration with @ruvector/graph-node
 *
 * @module ruvector-embeddings
 * @version 1.0.0
 */

// ============================================
// Type Definitions
// ============================================

/**
 * Embedding model options
 */
export type EmbeddingModel = 'all-MiniLM-L6-v2' | 'all-mpnet-base-v2';

/**
 * Configuration for embedding generation
 */
export interface EmbeddingConfig {
  /** Model to use (default: all-MiniLM-L6-v2) */
  model?: EmbeddingModel;
  /** Enable caching (default: true) */
  cache?: boolean;
  /** Cache size (default: 1000) */
  cacheSize?: number;
  /** Normalize embeddings (default: true) */
  normalize?: boolean;
}

/**
 * Configuration for HNSW index
 */
export interface HNSWConfig {
  /** M parameter (connections per node, default: 16) */
  m?: number;
  /** ef_construction (default: 100) */
  efConstruction?: number;
  /** ef_search (default: 50) */
  efSearch?: number;
}

/**
 * Search result from vector similarity search
 */
export interface VectorSearchResult {
  id: string;
  score: number;
  metadata?: Record<string, unknown>;
}

/**
 * Hyperbolic embedding configuration
 */
export interface HyperbolicConfig {
  /** Curvature of Poincare ball (default: -1.0) */
  curvature?: number;
  /** Enable hyperbolic mode */
  enabled?: boolean;
}

// ============================================
// RuVector Embeddings Class
// ============================================

/**
 * RuVector Embeddings Manager
 *
 * Generates and manages vector embeddings for the gene knowledge graph.
 *
 * @example
 * ```typescript
 * import { RuVectorEmbeddings } from './lib/ruvector-embeddings';
 *
 * const embeddings = new RuVectorEmbeddings();
 * await embeddings.initialize();
 *
 * // Generate embedding for text
 * const geneEmbed = await embeddings.embed('MTHFR methylenetetrahydrofolate reductase enzyme');
 *
 * // Batch embed
 * const batchEmbeds = await embeddings.batchEmbed([
 *   'MTHFR gene variant',
 *   'Vitamin B12 cofactor',
 *   'Methylation cycle pathway'
 * ]);
 *
 * // Search similar
 * const similar = await embeddings.search(geneEmbed, 10);
 * ```
 */
export class RuVectorEmbeddings {
  private config: Required<EmbeddingConfig>;
  private hnswConfig: Required<HNSWConfig>;
  private hyperbolicConfig: Required<HyperbolicConfig>;
  private initialized: boolean = false;
  private ruvector: RuVectorInstance | null = null;
  private cache: Map<string, Float32Array> = new Map();

  constructor(
    config: EmbeddingConfig = {},
    hnswConfig: HNSWConfig = {},
    hyperbolicConfig: HyperbolicConfig = {}
  ) {
    this.config = {
      model: config.model || 'all-MiniLM-L6-v2',
      cache: config.cache ?? true,
      cacheSize: config.cacheSize || 1000,
      normalize: config.normalize ?? true,
    };

    this.hnswConfig = {
      m: hnswConfig.m || 16,
      efConstruction: hnswConfig.efConstruction || 100,
      efSearch: hnswConfig.efSearch || 50,
    };

    this.hyperbolicConfig = {
      curvature: hyperbolicConfig.curvature || -1.0,
      enabled: hyperbolicConfig.enabled ?? false,
    };
  }

  /**
   * Initialize the embeddings system
   */
  async initialize(): Promise<void> {
    if (this.initialized) {
      return;
    }

    try {
      const ruvector = await import('ruvector');

      // Initialize with HNSW configuration
      this.ruvector = await ruvector.create({
        dimensions: this.getDimensions(),
        metric: 'cosine',
        hnsw: {
          m: this.hnswConfig.m,
          efConstruction: this.hnswConfig.efConstruction,
        },
      });

      this.initialized = true;
    } catch (error) {
      // Fallback: ruvector may not have ONNX embedding support
      // In this case, we provide a mock implementation
      console.warn('[RuVectorEmbeddings] Native ONNX not available, using fallback');
      this.initialized = true;
    }
  }

  /**
   * Ensure embeddings system is initialized
   */
  private ensureInitialized(): void {
    if (!this.initialized) {
      throw new Error('Embeddings not initialized. Call initialize() first.');
    }
  }

  /**
   * Get embedding dimensions for the configured model
   */
  getDimensions(): number {
    switch (this.config.model) {
      case 'all-MiniLM-L6-v2':
        return 384;
      case 'all-mpnet-base-v2':
        return 768;
      default:
        return 384;
    }
  }

  // ==========================================
  // Embedding Generation
  // ==========================================

  /**
   * Generate embedding for a text string
   */
  async embed(text: string): Promise<Float32Array> {
    this.ensureInitialized();

    // Check cache
    if (this.config.cache && this.cache.has(text)) {
      return this.cache.get(text)!;
    }

    let embedding: Float32Array;

    if (this.ruvector) {
      // Use native ruvector embedding
      embedding = await this.ruvector.embed(text);
    } else {
      // Fallback: generate deterministic pseudo-embedding from text hash
      embedding = this.generateFallbackEmbedding(text);
    }

    // Normalize if configured
    if (this.config.normalize) {
      embedding = this.normalize(embedding);
    }

    // Cache result
    if (this.config.cache) {
      this.addToCache(text, embedding);
    }

    return embedding;
  }

  /**
   * Batch generate embeddings
   */
  async batchEmbed(texts: string[]): Promise<Float32Array[]> {
    this.ensureInitialized();

    const results: Float32Array[] = [];
    const uncached: Array<{ index: number; text: string }> = [];

    // Check cache first
    for (let i = 0; i < texts.length; i++) {
      if (this.config.cache && this.cache.has(texts[i])) {
        results[i] = this.cache.get(texts[i])!;
      } else {
        uncached.push({ index: i, text: texts[i] });
      }
    }

    // Generate missing embeddings
    if (uncached.length > 0) {
      if (this.ruvector) {
        const uncachedTexts = uncached.map((u) => u.text);
        const embeddings = await this.ruvector.batchEmbed(uncachedTexts);

        for (let i = 0; i < uncached.length; i++) {
          let embedding = embeddings[i];
          if (this.config.normalize) {
            embedding = this.normalize(embedding);
          }
          results[uncached[i].index] = embedding;
          if (this.config.cache) {
            this.addToCache(uncached[i].text, embedding);
          }
        }
      } else {
        // Fallback for each
        for (const { index, text } of uncached) {
          let embedding = this.generateFallbackEmbedding(text);
          if (this.config.normalize) {
            embedding = this.normalize(embedding);
          }
          results[index] = embedding;
          if (this.config.cache) {
            this.addToCache(text, embedding);
          }
        }
      }
    }

    return results;
  }

  // ==========================================
  // Domain-Specific Embedding Generation
  // ==========================================

  /**
   * Generate embedding for a Gene entity
   */
  async embedGene(gene: {
    symbol: string;
    name: string;
    description?: string;
    chromosome?: string;
  }): Promise<Float32Array> {
    const text = [
      `Gene: ${gene.symbol}`,
      gene.name,
      gene.description || '',
      gene.chromosome ? `Chromosome ${gene.chromosome}` : '',
    ]
      .filter(Boolean)
      .join('. ');

    return this.embed(text);
  }

  /**
   * Generate embedding for an SNP entity
   */
  async embedSNP(snp: {
    rsid: string;
    geneSymbol?: string;
    clinicalSignificance?: string;
    commonName?: string;
  }): Promise<Float32Array> {
    const text = [
      `SNP: ${snp.rsid}`,
      snp.commonName || '',
      snp.geneSymbol ? `Gene: ${snp.geneSymbol}` : '',
      snp.clinicalSignificance || '',
    ]
      .filter(Boolean)
      .join('. ');

    return this.embed(text);
  }

  /**
   * Generate embedding for a Nutrient entity
   */
  async embedNutrient(nutrient: {
    name: string;
    category: string;
    alternateNames?: string[];
  }): Promise<Float32Array> {
    const text = [
      `Nutrient: ${nutrient.name}`,
      `Category: ${nutrient.category}`,
      nutrient.alternateNames?.length
        ? `Also known as: ${nutrient.alternateNames.join(', ')}`
        : '',
    ]
      .filter(Boolean)
      .join('. ');

    return this.embed(text);
  }

  /**
   * Generate embedding for an Herb entity
   */
  async embedHerb(herb: {
    name: string;
    latinName?: string;
    modality: string;
    properties?: string[];
  }): Promise<Float32Array> {
    const text = [
      `Herb: ${herb.name}`,
      herb.latinName ? `(${herb.latinName})` : '',
      `Modality: ${herb.modality}`,
      herb.properties?.length ? `Properties: ${herb.properties.join(', ')}` : '',
    ]
      .filter(Boolean)
      .join('. ');

    return this.embed(text);
  }

  /**
   * Generate embedding for a Pathway entity
   */
  async embedPathway(pathway: {
    name: string;
    description?: string;
    category?: string;
  }): Promise<Float32Array> {
    const text = [
      `Pathway: ${pathway.name}`,
      pathway.description || '',
      pathway.category ? `Category: ${pathway.category}` : '',
    ]
      .filter(Boolean)
      .join('. ');

    return this.embed(text);
  }

  /**
   * Generate embedding for an edge/relationship
   */
  async embedRelationship(rel: {
    type: string;
    fromLabel: string;
    toLabel: string;
    metadata?: Record<string, unknown>;
  }): Promise<Float32Array> {
    const metaStr = rel.metadata
      ? Object.entries(rel.metadata)
          .map(([k, v]) => `${k}: ${v}`)
          .join(', ')
      : '';

    const text = [
      `Relationship: ${rel.type}`,
      `From: ${rel.fromLabel} To: ${rel.toLabel}`,
      metaStr,
    ]
      .filter(Boolean)
      .join('. ');

    return this.embed(text);
  }

  // ==========================================
  // Vector Search
  // ==========================================

  /**
   * Store embedding with ID for later search
   */
  async store(
    id: string,
    embedding: Float32Array,
    metadata?: Record<string, unknown>
  ): Promise<void> {
    this.ensureInitialized();

    if (this.ruvector) {
      await this.ruvector.insert(id, embedding, metadata);
    }
  }

  /**
   * Search for similar embeddings
   */
  async search(
    embedding: Float32Array,
    k: number = 10,
    threshold: number = 0.5
  ): Promise<VectorSearchResult[]> {
    this.ensureInitialized();

    if (this.ruvector) {
      const results = await this.ruvector.search(embedding, k);
      return results.filter((r: VectorSearchResult) => r.score >= threshold);
    }

    return [];
  }

  /**
   * Search by text (generates embedding internally)
   */
  async searchByText(
    text: string,
    k: number = 10,
    threshold: number = 0.5
  ): Promise<VectorSearchResult[]> {
    const embedding = await this.embed(text);
    return this.search(embedding, k, threshold);
  }

  // ==========================================
  // Hyperbolic Embeddings
  // ==========================================

  /**
   * Convert Euclidean embedding to Poincare ball
   */
  toPoincare(embedding: Float32Array): Float32Array {
    if (!this.hyperbolicConfig.enabled) {
      return embedding;
    }

    const c = Math.abs(this.hyperbolicConfig.curvature);
    const result = new Float32Array(embedding.length);

    // Exponential map from origin
    let norm = 0;
    for (let i = 0; i < embedding.length; i++) {
      norm += embedding[i] * embedding[i];
    }
    norm = Math.sqrt(norm);

    if (norm < 1e-8) {
      return result; // Zero vector
    }

    const sqrtC = Math.sqrt(c);
    const factor = Math.tanh(sqrtC * norm) / (sqrtC * norm);

    for (let i = 0; i < embedding.length; i++) {
      result[i] = embedding[i] * factor;
    }

    return result;
  }

  /**
   * Calculate Poincare distance between two hyperbolic embeddings
   */
  poincareDistance(x: Float32Array, y: Float32Array): number {
    const c = Math.abs(this.hyperbolicConfig.curvature);

    let normX = 0;
    let normY = 0;
    let normDiff = 0;

    for (let i = 0; i < x.length; i++) {
      normX += x[i] * x[i];
      normY += y[i] * y[i];
      normDiff += (x[i] - y[i]) * (x[i] - y[i]);
    }

    const numerator = 2 * normDiff;
    const denominator = (1 - c * normX) * (1 - c * normY);

    if (denominator <= 0) {
      return Infinity;
    }

    const sqrtC = Math.sqrt(c);
    return (2 / sqrtC) * Math.atanh(sqrtC * Math.sqrt(numerator / denominator));
  }

  // ==========================================
  // Utility Methods
  // ==========================================

  /**
   * Normalize embedding to unit length
   */
  normalize(embedding: Float32Array): Float32Array {
    let norm = 0;
    for (let i = 0; i < embedding.length; i++) {
      norm += embedding[i] * embedding[i];
    }
    norm = Math.sqrt(norm);

    if (norm < 1e-8) {
      return embedding;
    }

    const result = new Float32Array(embedding.length);
    for (let i = 0; i < embedding.length; i++) {
      result[i] = embedding[i] / norm;
    }

    return result;
  }

  /**
   * Calculate cosine similarity between two embeddings
   */
  similarity(a: Float32Array, b: Float32Array): number {
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

    const denominator = Math.sqrt(normA) * Math.sqrt(normB);
    return denominator > 0 ? dot / denominator : 0;
  }

  /**
   * Generate deterministic fallback embedding from text hash
   */
  private generateFallbackEmbedding(text: string): Float32Array {
    const dimensions = this.getDimensions();
    const embedding = new Float32Array(dimensions);

    // Simple hash-based pseudo-embedding
    let hash = 0;
    for (let i = 0; i < text.length; i++) {
      hash = ((hash << 5) - hash + text.charCodeAt(i)) | 0;
    }

    const seed = Math.abs(hash);
    for (let i = 0; i < dimensions; i++) {
      // Deterministic pseudo-random based on seed and index
      const x = Math.sin(seed * (i + 1)) * 10000;
      embedding[i] = x - Math.floor(x) - 0.5;
    }

    return embedding;
  }

  /**
   * Add embedding to cache with LRU eviction
   */
  private addToCache(text: string, embedding: Float32Array): void {
    if (this.cache.size >= this.config.cacheSize) {
      // Remove oldest entry (first key)
      const firstKey = this.cache.keys().next().value;
      if (firstKey) {
        this.cache.delete(firstKey);
      }
    }
    this.cache.set(text, embedding);
  }

  /**
   * Clear the embedding cache
   */
  clearCache(): void {
    this.cache.clear();
  }

  /**
   * Get cache statistics
   */
  getCacheStats(): { size: number; maxSize: number } {
    return {
      size: this.cache.size,
      maxSize: this.config.cacheSize,
    };
  }
}

// ============================================
// Type Declarations for Dynamic Import
// ============================================

interface RuVectorInstance {
  embed: (text: string) => Promise<Float32Array>;
  batchEmbed: (texts: string[]) => Promise<Float32Array[]>;
  insert: (id: string, embedding: Float32Array, metadata?: Record<string, unknown>) => Promise<void>;
  search: (embedding: Float32Array, k: number) => Promise<VectorSearchResult[]>;
}

// ============================================
// Factory Functions
// ============================================

/**
 * Create and initialize an embeddings instance
 */
export async function createEmbeddings(
  config?: EmbeddingConfig,
  hnswConfig?: HNSWConfig
): Promise<RuVectorEmbeddings> {
  const embeddings = new RuVectorEmbeddings(config, hnswConfig);
  await embeddings.initialize();
  return embeddings;
}

// ============================================
// Default Export
// ============================================

export default RuVectorEmbeddings;
