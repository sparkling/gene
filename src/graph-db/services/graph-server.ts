/**
 * RuVector Graph Database Server
 *
 * REST API server for THREE WORLDS knowledge graph
 * Using @ruvector/graph-node native bindings
 *
 * THREE WORLDS Architecture:
 * - World 1: Genetics (Gene, SNP, Pathway)
 * - World 2: Traditional Medicine (Herb, Formula, Modality)
 * - World 3: Nutrition (Nutrient, Food, Biomarker)
 */

import express, { Request, Response, NextFunction } from 'express';
import cors from 'cors';
import compression from 'compression';
import helmet from 'helmet';
import rateLimit from 'express-rate-limit';
import { v4 as uuidv4 } from 'uuid';

// Import RuVector Graph Database
const { GraphDatabase, JsDistanceMetric } = require('@ruvector/graph-node');

// ============================================
// Configuration
// ============================================

const config = {
  port: parseInt(process.env.PORT || '3100', 10),
  storagePath: process.env.GRAPH_STORAGE_PATH || './data/knowledge.db',
  dimensions: parseInt(process.env.GRAPH_DIMENSIONS || '384', 10),
  distanceMetric: (process.env.GRAPH_DISTANCE_METRIC || 'Cosine') as 'Cosine' | 'Euclidean' | 'DotProduct',
  logLevel: process.env.LOG_LEVEL || 'info',
};

// ============================================
// Graph Database Instance
// ============================================

let db: InstanceType<typeof GraphDatabase>;

function initializeDatabase(): void {
  console.log(`[GraphServer] Initializing database at: ${config.storagePath}`);
  console.log(`[GraphServer] Dimensions: ${config.dimensions}, Metric: ${config.distanceMetric}`);

  db = new GraphDatabase({
    distanceMetric: config.distanceMetric,
    dimensions: config.dimensions,
    storagePath: config.storagePath,
  });

  console.log(`[GraphServer] Database initialized. Persistent: ${db.isPersistent()}`);
}

// ============================================
// THREE WORLDS Node Types
// ============================================

type World = 'genetics' | 'traditional_medicine' | 'nutrition';

interface NodeLabels {
  genetics: ['Gene', 'SNP', 'Pathway', 'Condition', 'Publication'];
  traditional_medicine: ['Herb', 'Formula', 'Modality', 'TCMPattern', 'Meridian'];
  nutrition: ['Nutrient', 'Food', 'Biomarker', 'Supplement', 'FoodGroup'];
}

const WORLD_LABELS: NodeLabels = {
  genetics: ['Gene', 'SNP', 'Pathway', 'Condition', 'Publication'],
  traditional_medicine: ['Herb', 'Formula', 'Modality', 'TCMPattern', 'Meridian'],
  nutrition: ['Nutrient', 'Food', 'Biomarker', 'Supplement', 'FoodGroup'],
};

// Cross-world relationship types
const CROSS_WORLD_EDGES = [
  'AFFECTS',           // SNP -> Nutrient (genetics -> nutrition)
  'TREATS',            // Herb -> Condition (traditional -> genetics)
  'CONTAINS',          // Food -> Nutrient (nutrition internal)
  'HAS_VARIANT',       // Gene -> SNP (genetics internal)
  'INVOLVED_IN',       // Gene -> Pathway (genetics internal)
  'INGREDIENT_OF',     // Herb -> Formula (traditional internal)
  'REQUIRES_COFACTOR', // Pathway -> Nutrient (genetics -> nutrition)
  'ASSOCIATED_WITH',   // SNP -> Condition (genetics internal)
  'HAS_RESEARCH',      // Any -> Publication
  'TARGETS',           // Herb -> Pathway (traditional -> genetics)
  'MODULATES',         // Nutrient -> Gene (nutrition -> genetics)
  'DEFICIENCY_CAUSES', // Nutrient -> Condition (nutrition -> genetics)
];

// ============================================
// Express Application
// ============================================

const app = express();

// Middleware
app.use(helmet());
app.use(cors());
app.use(compression());
app.use(express.json({ limit: '10mb' }));

// Rate limiting
const limiter = rateLimit({
  windowMs: 60 * 1000, // 1 minute
  max: 1000, // 1000 requests per minute
  standardHeaders: true,
  legacyHeaders: false,
});
app.use(limiter);

// Request logging
app.use((req: Request, _res: Response, next: NextFunction) => {
  if (config.logLevel === 'debug') {
    console.log(`[${new Date().toISOString()}] ${req.method} ${req.path}`);
  }
  next();
});

// ============================================
// Health & Status Endpoints
// ============================================

app.get('/health', async (_req: Request, res: Response) => {
  try {
    const stats = await db.stats();
    res.json({
      status: 'healthy',
      persistent: db.isPersistent(),
      storagePath: db.getStoragePath(),
      stats,
      timestamp: new Date().toISOString(),
    });
  } catch (error) {
    res.status(503).json({ status: 'unhealthy', error: String(error) });
  }
});

app.get('/stats', async (_req: Request, res: Response) => {
  try {
    const stats = await db.stats();
    res.json({
      ...stats,
      worlds: WORLD_LABELS,
      crossWorldEdges: CROSS_WORLD_EDGES,
    });
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Node Operations
// ============================================

/**
 * Create a node in the graph
 * POST /nodes
 * Body: { id, embedding, labels, properties }
 */
app.post('/nodes', async (req: Request, res: Response) => {
  try {
    const { id, embedding, labels, properties } = req.body;

    if (!id || !embedding) {
      return res.status(400).json({ error: 'id and embedding are required' });
    }

    // Convert embedding to Float32Array if needed
    const embeddingArray = embedding instanceof Float32Array
      ? embedding
      : new Float32Array(embedding);

    if (embeddingArray.length !== config.dimensions) {
      return res.status(400).json({
        error: `Embedding must have ${config.dimensions} dimensions, got ${embeddingArray.length}`,
      });
    }

    const nodeId = await db.createNode({
      id,
      embedding: embeddingArray,
      labels: labels || [],
      properties: properties || {},
    });

    res.status(201).json({ nodeId, success: true });
  } catch (error) {
    console.error('[CreateNode Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Batch create nodes
 * POST /nodes/batch
 * Body: { nodes: [{ id, embedding, labels, properties }] }
 */
app.post('/nodes/batch', async (req: Request, res: Response) => {
  try {
    const { nodes } = req.body;

    if (!Array.isArray(nodes) || nodes.length === 0) {
      return res.status(400).json({ error: 'nodes array is required' });
    }

    const preparedNodes = nodes.map((node: any) => ({
      id: node.id,
      embedding: node.embedding instanceof Float32Array
        ? node.embedding
        : new Float32Array(node.embedding),
      labels: node.labels || [],
      properties: node.properties || {},
    }));

    const result = await db.batchInsert({
      nodes: preparedNodes,
      edges: [],
    });

    res.status(201).json({
      success: true,
      nodeIds: result.nodeIds,
      count: result.nodeIds.length,
    });
  } catch (error) {
    console.error('[BatchCreateNodes Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Edge Operations
// ============================================

/**
 * Create an edge between nodes
 * POST /edges
 * Body: { from, to, description, embedding, confidence, metadata }
 */
app.post('/edges', async (req: Request, res: Response) => {
  try {
    const { from, to, description, embedding, confidence, metadata } = req.body;

    if (!from || !to || !description || !embedding) {
      return res.status(400).json({
        error: 'from, to, description, and embedding are required',
      });
    }

    const embeddingArray = embedding instanceof Float32Array
      ? embedding
      : new Float32Array(embedding);

    const edgeId = await db.createEdge({
      from,
      to,
      description,
      embedding: embeddingArray,
      confidence: confidence ?? 1.0,
      metadata: metadata || {},
    });

    res.status(201).json({ edgeId, success: true });
  } catch (error) {
    console.error('[CreateEdge Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Batch create edges
 * POST /edges/batch
 */
app.post('/edges/batch', async (req: Request, res: Response) => {
  try {
    const { edges } = req.body;

    if (!Array.isArray(edges) || edges.length === 0) {
      return res.status(400).json({ error: 'edges array is required' });
    }

    const preparedEdges = edges.map((edge: any) => ({
      from: edge.from,
      to: edge.to,
      description: edge.description,
      embedding: edge.embedding instanceof Float32Array
        ? edge.embedding
        : new Float32Array(edge.embedding),
      confidence: edge.confidence ?? 1.0,
      metadata: edge.metadata || {},
    }));

    const result = await db.batchInsert({
      nodes: [],
      edges: preparedEdges,
    });

    res.status(201).json({
      success: true,
      edgeIds: result.edgeIds,
      count: result.edgeIds.length,
    });
  } catch (error) {
    console.error('[BatchCreateEdges Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Hyperedge Operations (Multi-node relationships)
// ============================================

/**
 * Create a hyperedge connecting multiple nodes
 * POST /hyperedges
 * Body: { nodes, description, embedding, confidence, metadata }
 *
 * Use case: Formula compositions, pathway complexes, gene-nutrient-condition triangles
 */
app.post('/hyperedges', async (req: Request, res: Response) => {
  try {
    const { nodes, description, embedding, confidence, metadata } = req.body;

    if (!Array.isArray(nodes) || nodes.length < 2 || !description || !embedding) {
      return res.status(400).json({
        error: 'nodes (array of 2+), description, and embedding are required',
      });
    }

    const embeddingArray = embedding instanceof Float32Array
      ? embedding
      : new Float32Array(embedding);

    const hyperedgeId = await db.createHyperedge({
      nodes,
      description,
      embedding: embeddingArray,
      confidence: confidence ?? 1.0,
      metadata: metadata || {},
    });

    res.status(201).json({ hyperedgeId, success: true });
  } catch (error) {
    console.error('[CreateHyperedge Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Search hyperedges by similarity
 * POST /hyperedges/search
 * Body: { embedding, k }
 */
app.post('/hyperedges/search', async (req: Request, res: Response) => {
  try {
    const { embedding, k = 10 } = req.body;

    if (!embedding) {
      return res.status(400).json({ error: 'embedding is required' });
    }

    const embeddingArray = embedding instanceof Float32Array
      ? embedding
      : new Float32Array(embedding);

    const results = await db.searchHyperedges({
      embedding: embeddingArray,
      k,
    });

    res.json({ results, count: results.length });
  } catch (error) {
    console.error('[SearchHyperedges Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Query Operations (Cypher-like)
// ============================================

/**
 * Execute a Cypher query
 * POST /query
 * Body: { cypher }
 */
app.post('/query', async (req: Request, res: Response) => {
  try {
    const { cypher } = req.body;

    if (!cypher) {
      return res.status(400).json({ error: 'cypher query is required' });
    }

    const result = await db.query(cypher);
    res.json(result);
  } catch (error) {
    console.error('[Query Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Execute a Cypher query (synchronous)
 * POST /query/sync
 * Body: { cypher }
 */
app.post('/query/sync', (req: Request, res: Response) => {
  try {
    const { cypher } = req.body;

    if (!cypher) {
      return res.status(400).json({ error: 'cypher query is required' });
    }

    const result = db.querySync(cypher);
    res.json(result);
  } catch (error) {
    console.error('[QuerySync Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Graph Traversal Operations
// ============================================

/**
 * Get k-hop neighbors from a starting node
 * GET /neighbors/:nodeId?k=2
 */
app.get('/neighbors/:nodeId', async (req: Request, res: Response) => {
  try {
    const { nodeId } = req.params;
    const k = parseInt(req.query.k as string || '1', 10);

    const neighbors = await db.kHopNeighbors(nodeId, k);
    res.json({ nodeId, k, neighbors, count: neighbors.length });
  } catch (error) {
    console.error('[KHopNeighbors Error]', error);
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// THREE WORLDS Specific Queries
// ============================================

/**
 * Find SNPs affecting a nutrient
 * GET /worlds/snp-nutrient/:nutrientId
 */
app.get('/worlds/snp-nutrient/:nutrientId', async (req: Request, res: Response) => {
  try {
    const { nutrientId } = req.params;
    const cypher = `
      MATCH (s:SNP)-[a:AFFECTS]->(n:Nutrient)
      WHERE n.id = '${nutrientId}'
      RETURN s, a
    `;
    const result = await db.query(cypher);
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Find herbs treating a condition via gene pathways
 * GET /worlds/herb-condition/:conditionId
 */
app.get('/worlds/herb-condition/:conditionId', async (req: Request, res: Response) => {
  try {
    const { conditionId } = req.params;
    const cypher = `
      MATCH (h:Herb)-[:TARGETS]->(p:Pathway)<-[:INVOLVED_IN]-(g:Gene)-[:HAS_VARIANT]->(s:SNP)-[:ASSOCIATED_WITH]->(c:Condition)
      WHERE c.id = '${conditionId}'
      RETURN h, p, g, s, c
    `;
    const result = await db.query(cypher);
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Find nutrient-gene-condition pathways
 * GET /worlds/nutrient-gene/:geneSymbol
 */
app.get('/worlds/nutrient-gene/:geneSymbol', async (req: Request, res: Response) => {
  try {
    const { geneSymbol } = req.params;
    const cypher = `
      MATCH (n:Nutrient)-[:MODULATES]->(g:Gene)
      WHERE g.symbol = '${geneSymbol}'
      RETURN n, g
    `;
    const result = await db.query(cypher);
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Get formula composition (hyperedge)
 * GET /worlds/formula/:formulaId
 */
app.get('/worlds/formula/:formulaId', async (req: Request, res: Response) => {
  try {
    const { formulaId } = req.params;
    const cypher = `
      MATCH (h:Herb)-[i:INGREDIENT_OF]->(f:Formula)
      WHERE f.id = '${formulaId}'
      RETURN h, i, f
    `;
    const result = await db.query(cypher);
    res.json(result);
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Transaction Operations
// ============================================

/**
 * Begin a transaction
 * POST /transactions/begin
 */
app.post('/transactions/begin', async (_req: Request, res: Response) => {
  try {
    const txId = await db.begin();
    res.json({ txId, success: true });
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Commit a transaction
 * POST /transactions/:txId/commit
 */
app.post('/transactions/:txId/commit', async (req: Request, res: Response) => {
  try {
    const { txId } = req.params;
    await db.commit(txId);
    res.json({ success: true });
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

/**
 * Rollback a transaction
 * POST /transactions/:txId/rollback
 */
app.post('/transactions/:txId/rollback', async (req: Request, res: Response) => {
  try {
    const { txId } = req.params;
    await db.rollback(txId);
    res.json({ success: true });
  } catch (error) {
    res.status(500).json({ error: String(error) });
  }
});

// ============================================
// Error Handling
// ============================================

app.use((err: Error, _req: Request, res: Response, _next: NextFunction) => {
  console.error('[Server Error]', err);
  res.status(500).json({
    error: 'Internal server error',
    message: err.message,
  });
});

// 404 handler
app.use((_req: Request, res: Response) => {
  res.status(404).json({ error: 'Not found' });
});

// ============================================
// Server Startup
// ============================================

async function startServer(): Promise<void> {
  try {
    // Initialize database
    initializeDatabase();

    // Start server
    app.listen(config.port, () => {
      console.log(`[GraphServer] Server running on port ${config.port}`);
      console.log(`[GraphServer] Health check: http://localhost:${config.port}/health`);
      console.log(`[GraphServer] Storage: ${config.storagePath}`);
    });

    // Graceful shutdown
    process.on('SIGTERM', async () => {
      console.log('[GraphServer] Shutting down...');
      process.exit(0);
    });

    process.on('SIGINT', async () => {
      console.log('[GraphServer] Interrupted, shutting down...');
      process.exit(0);
    });
  } catch (error) {
    console.error('[GraphServer] Failed to start:', error);
    process.exit(1);
  }
}

startServer();

export { app, db };
