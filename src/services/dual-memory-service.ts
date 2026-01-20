/**
 * Dual Memory Service
 *
 * Two separate RuVector databases:
 * - Operational: Claude Flow internals (expendable)
 * - User: Your graph database (critical)
 */

import {
  STORAGE_CONFIG,
  OperationalNamespace,
  UserNamespace,
  DatabaseConfig
} from '../config/storage-config';

export interface StoreOptions {
  type?: string;
  tags?: string[];
  metadata?: Record<string, any>;
}

export interface SearchOptions {
  namespace?: string;
  limit?: number;
  threshold?: number;
}

export interface SearchResult {
  id: string;
  score: number;
  value: any;
  metadata: Record<string, any>;
}

/**
 * RuVector Database Wrapper
 */
class RuVectorStore {
  private db: any = null;
  private config: DatabaseConfig;
  private namespaces: Record<string, string>;
  private initialized = false;

  constructor(config: DatabaseConfig, namespaces: Record<string, string>) {
    this.config = config;
    this.namespaces = namespaces;
  }

  async initialize(): Promise<void> {
    if (this.initialized) return;

    try {
      const ruvector = await import('ruvector').catch(() => null);

      if (ruvector?.VectorDB) {
        this.db = new ruvector.VectorDB({
          dimensions: this.config.dimensions,
          storagePath: this.config.storagePath,
          distanceMetric: 'cosine',
          hnswConfig: {
            m: this.config.hnsw.M,
            efConstruction: this.config.hnsw.efConstruction,
            efSearch: this.config.hnsw.efSearch,
          },
        });
        console.log(`[RuVectorStore] Initialized: ${this.config.storagePath}`);
      }

      this.initialized = true;
    } catch (error) {
      console.warn('[RuVectorStore] Fallback mode:', error);
      this.initialized = true;
    }
  }

  private getPrefix(namespace: string): string {
    return this.namespaces[namespace] || 'default:';
  }

  private generateId(): string {
    return `${Date.now()}_${Math.random().toString(36).slice(2, 8)}`;
  }

  async store(namespace: string, key: string, value: any, embedding: Float32Array, options?: StoreOptions): Promise<string> {
    await this.initialize();

    const prefix = this.getPrefix(namespace);
    const id = `${prefix}${key || this.generateId()}`;

    if (this.db) {
      await this.db.insert({
        id,
        vector: embedding,
        metadata: {
          value,
          type: options?.type || 'default',
          tags: options?.tags || [],
          ...options?.metadata,
          _namespace: namespace,
          _created: Date.now(),
        }
      });
    }

    return id;
  }

  async search(query: Float32Array, namespace?: string, limit = 10): Promise<SearchResult[]> {
    await this.initialize();

    if (!this.db) return [];

    const results = await this.db.search({
      vector: query,
      k: limit * 2,  // Over-fetch to filter
      efSearch: this.config.hnsw.efSearch,
    });

    let filtered = results;

    if (namespace) {
      const prefix = this.getPrefix(namespace);
      filtered = results.filter((r: any) => r.id.startsWith(prefix));
    }

    return filtered.slice(0, limit).map((r: any) => ({
      id: r.id,
      score: r.score,
      value: r.metadata?.value,
      metadata: r.metadata,
    }));
  }

  async get(id: string): Promise<any | null> {
    await this.initialize();
    if (!this.db) return null;

    const result = await this.db.get(id);
    return result?.metadata?.value || null;
  }

  async delete(id: string): Promise<boolean> {
    await this.initialize();
    if (!this.db) return false;

    await this.db.delete(id);
    return true;
  }

  async save(): Promise<void> {
    if (this.db?.save) {
      await this.db.save();
    }
  }

  async count(): Promise<number> {
    await this.initialize();
    if (!this.db) return 0;
    return this.db.count?.() || 0;
  }
}

/**
 * Dual Memory Service
 * Manages both operational and user databases
 */
export class DualMemoryService {
  private operational: RuVectorStore;
  private user: RuVectorStore;
  private embedder: any = null;

  constructor() {
    this.operational = new RuVectorStore(
      STORAGE_CONFIG.operational.db,
      STORAGE_CONFIG.operational.namespaces
    );
    this.user = new RuVectorStore(
      STORAGE_CONFIG.user.db,
      STORAGE_CONFIG.user.namespaces
    );
  }

  private async getEmbedder(): Promise<any> {
    if (this.embedder) return this.embedder;

    try {
      const ruvector = await import('ruvector').catch(() => null);
      if (ruvector?.IntelligenceEngine) {
        this.embedder = new ruvector.IntelligenceEngine({
          embeddingDim: 768,
          enableOnnx: true,
        });
      }
    } catch (e) {
      console.warn('[DualMemoryService] Embedder not available');
    }

    return this.embedder;
  }

  private async embed(text: string, dimensions: number): Promise<Float32Array> {
    const embedder = await this.getEmbedder();

    if (embedder) {
      const embedding = await embedder.embedAsync(text);
      // Resize if needed
      if (embedding.length !== dimensions) {
        const resized = new Float32Array(dimensions);
        for (let i = 0; i < Math.min(embedding.length, dimensions); i++) {
          resized[i] = embedding[i];
        }
        return resized;
      }
      return new Float32Array(embedding);
    }

    // Fallback: random embedding (for testing)
    const fallback = new Float32Array(dimensions);
    for (let i = 0; i < dimensions; i++) {
      fallback[i] = Math.random() - 0.5;
    }
    return fallback;
  }

  // =========================================================================
  // OPERATIONAL DATABASE (Claude Flow internals)
  // =========================================================================

  async storeOperational(
    namespace: OperationalNamespace,
    key: string,
    value: any,
    options?: StoreOptions
  ): Promise<string> {
    const text = typeof value === 'string' ? value : JSON.stringify(value);
    const embedding = await this.embed(text, STORAGE_CONFIG.operational.db.dimensions);
    return this.operational.store(namespace, key, value, embedding, options);
  }

  async searchOperational(
    query: string,
    namespace?: OperationalNamespace,
    limit = 10
  ): Promise<SearchResult[]> {
    const embedding = await this.embed(query, STORAGE_CONFIG.operational.db.dimensions);
    return this.operational.search(embedding, namespace, limit);
  }

  async getOperational(id: string): Promise<any | null> {
    return this.operational.get(id);
  }

  async deleteOperational(id: string): Promise<boolean> {
    return this.operational.delete(id);
  }

  // =========================================================================
  // USER DATABASE (Your graph data - CRITICAL)
  // =========================================================================

  async storeUser(
    namespace: UserNamespace,
    key: string,
    value: any,
    options?: StoreOptions
  ): Promise<string> {
    const text = typeof value === 'string' ? value : JSON.stringify(value);
    const embedding = await this.embed(text, STORAGE_CONFIG.user.db.dimensions);
    return this.user.store(namespace, key, value, embedding, options);
  }

  async searchUser(
    query: string,
    namespace?: UserNamespace,
    limit = 10
  ): Promise<SearchResult[]> {
    const embedding = await this.embed(query, STORAGE_CONFIG.user.db.dimensions);
    return this.user.search(embedding, namespace, limit);
  }

  async getUser(id: string): Promise<any | null> {
    return this.user.get(id);
  }

  async deleteUser(id: string): Promise<boolean> {
    return this.user.delete(id);
  }

  // =========================================================================
  // CONVENIENCE METHODS
  // =========================================================================

  /** Store a document in user database */
  async storeDocument(key: string, content: any, metadata?: Record<string, any>): Promise<string> {
    return this.storeUser('documents', key, content, { metadata });
  }

  /** Store a graph entity */
  async storeEntity(key: string, entity: any, metadata?: Record<string, any>): Promise<string> {
    return this.storeUser('entities', key, entity, { type: 'entity', metadata });
  }

  /** Store a graph relationship */
  async storeRelationship(key: string, relationship: any, metadata?: Record<string, any>): Promise<string> {
    return this.storeUser('relationships', key, relationship, { type: 'relationship', metadata });
  }

  /** Log project activity (operational) */
  async logActivity(key: string, activity: any): Promise<string> {
    return this.storeOperational('project-activity', key, activity);
  }

  /** Store a learned pattern (operational) */
  async storePattern(key: string, pattern: any): Promise<string> {
    return this.storeOperational('learning-patterns', key, pattern);
  }

  // =========================================================================
  // PERSISTENCE
  // =========================================================================

  async saveAll(): Promise<void> {
    await Promise.all([
      this.operational.save(),
      this.user.save(),
    ]);
  }

  async getStats(): Promise<{ operational: number; user: number }> {
    return {
      operational: await this.operational.count(),
      user: await this.user.count(),
    };
  }
}

// Singleton
let instance: DualMemoryService | null = null;

export function getDualMemoryService(): DualMemoryService {
  if (!instance) {
    instance = new DualMemoryService();
  }
  return instance;
}

export default DualMemoryService;
