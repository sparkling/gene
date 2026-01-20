/**
 * Storage Configuration
 * Two separate RuVector databases:
 * 1. Operational - Claude Flow internals (expendable)
 * 2. User Data - Your graph database (critical)
 */

export interface HNSWConfig {
  M: number;
  efConstruction: number;
  efSearch: number;
}

export interface PersistenceConfig {
  autoPersist: boolean;
  persistInterval: number;
  walMode: boolean;
}

export interface DatabaseConfig {
  storagePath: string;
  dimensions: number;
  hnsw: HNSWConfig;
  quantization: 'none' | 'scalar' | 'binary' | 'product';
  persistence: PersistenceConfig;
}

// =============================================================================
// OPERATIONAL DATABASE - Claude Flow internals (expendable, can rebuild)
// =============================================================================
export const OPERATIONAL_DB: DatabaseConfig = {
  storagePath: '/home/claude/src/gene/.claude-flow/operational.db',
  dimensions: 384,  // Smaller for faster operations

  hnsw: {
    M: 16,              // Lower - operational data is small
    efConstruction: 200,
    efSearch: 100,
  },

  quantization: 'none',  // Small dataset, no need

  persistence: {
    autoPersist: true,
    persistInterval: 30000,  // 30 seconds - less critical
    walMode: false,
  }
};

export const OPERATIONAL_NAMESPACES = {
  'project-activity': 'proj:',   // What Claude does
  'learning-patterns': 'learn:', // SONA patterns
  'session-state': 'sess:',      // Conversation context
  'embeddings-cache': 'emb:',    // Cached embeddings
  'agent-routing': 'route:',     // Agent routing decisions
} as const;

// =============================================================================
// USER DATABASE - Your graph data (CRITICAL - cannot lose)
// =============================================================================
export const USER_DB: DatabaseConfig = {
  storagePath: '/home/claude/src/gene/.vectordb/user-data.db',
  dimensions: 768,  // Higher quality for user data

  hnsw: {
    M: 48,              // High - user data needs best recall
    efConstruction: 400,
    efSearch: 200,
  },

  quantization: 'scalar',  // 4x compression, 98-99% accuracy

  persistence: {
    autoPersist: true,
    persistInterval: 10000,  // 10 seconds - critical data
    walMode: true,           // Crash-safe
  }
};

export const USER_NAMESPACES = {
  'documents': 'doc:',       // Source documents
  'entities': 'ent:',        // Graph entities/nodes
  'relationships': 'rel:',   // Graph relationships/edges
  'annotations': 'ann:',     // User annotations
  'queries': 'qry:',         // Saved queries/searches
} as const;

// =============================================================================
// EXPORTS
// =============================================================================
export const STORAGE_CONFIG = {
  operational: {
    db: OPERATIONAL_DB,
    namespaces: OPERATIONAL_NAMESPACES,
  },
  user: {
    db: USER_DB,
    namespaces: USER_NAMESPACES,
  }
};

export type OperationalNamespace = keyof typeof OPERATIONAL_NAMESPACES;
export type UserNamespace = keyof typeof USER_NAMESPACES;
