# RuVector Unified Database Architecture Plan

**Document ID:** PLAN-RUVECTOR
**Status:** Under Review - Critical Findings
**Owner:** Claude
**Last Updated:** January 2026
**Version:** 3.0

---

## TL;DR

> ⚠️ **CRITICAL UPDATE (v3.0):** Hive-mind analysis of 7 specialist agents recommends **proceeding with caution** on RuVector. The Critical Analyst scored RuVector at **2.30/5.00** vs TypeDB/Apache AGE at **3.95/5.00** due to version immaturity (0.1.x), single maintainer risk, and critical open bugs.

**Revised Recommendation:**
- **Primary Option:** Evaluate **CozoDB** (Datalog + HNSW, Rust-native, MPL 2.0) or **Apache AGE** (PostgreSQL + Cypher)
- **If RuVector chosen:** Mandatory 3-week POC with escape plan documented
- **For biomedical SPARQL:** Consider **QLever** hybrid architecture

**Three-Database Architecture (unchanged):**
1. **Vector Database** - Operational data + vector embeddings (Claude Flow)
2. **Graph Database** - Knowledge graph for genes/SNPs/pathways/herbs relationships
3. **PostgreSQL (Separate)** - User domain data (HIPAA-compliant PII)

**Research Documentation:** See `docs/research/graphs/` (4,793 lines, 172 KB)

---

## Architecture Decision: RuVector Graph vs Neo4j

### Why RuVector Graph?

| Criteria | RuVector Graph | Neo4j |
|----------|----------------|-------|
| **Performance** | 131K ops/sec (native Rust) | ~50K ops/sec |
| **Vector-Graph Hybrid** | Native support | Requires plugin |
| **Cypher Queries** | Full compatibility | Native |
| **Hypergraphs** | Built-in | Not supported |
| **ACID Transactions** | Yes | Yes |
| **Stack Consistency** | Unified RuVector | Separate system |
| **Operational Cost** | Single stack | Additional infra |
| **Embedding Integration** | Native (same vectors) | External sync |

### Key RuVector Graph Features

From `@ruvector/graph-node`:

```typescript
// Hypergraph support for complex relationships
interface Hyperedge {
  id?: string;
  nodes: string[];  // Connect multiple nodes
  type: string;
  properties?: Record<string, any>;
}

// Full Cypher query support
cypher(query: string, params?: Record<string, any>): CypherResult;

// Graph algorithms built-in
pageRank(iterations?: number, dampingFactor?: number): Map<string, number>;
connectedComponents(): string[][];
communities(): Map<string, number>;  // Louvain algorithm
betweennessCentrality(): Map<string, number>;
shortestPath(from: string, to: string, maxDepth?: number): PathResult | null;
allPaths(from: string, to: string, maxDepth?: number, maxPaths?: number): PathResult[];
```

---

## Revised Three-Database Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         APPLICATION LAYER (Next.js/Node.js)                  │
└───────────────────────────────────┬─────────────────────────────────────────┘
                                    │
        ┌───────────────────────────┼───────────────────────────┐
        ▼                           ▼                           ▼
┌───────────────────┐   ┌────────────────────────┐   ┌────────────────────────┐
│  RuVector         │   │  RuVector Graph        │   │  PostgreSQL            │
│  PostgreSQL       │   │  (@ruvector/graph-node)│   │  (User Domain)         │
├───────────────────┤   ├────────────────────────┤   ├────────────────────────┤
│                   │   │                        │   │                        │
│  • Embeddings     │   │  • Gene nodes          │   │  • users               │
│  • Patterns       │   │  • SNP nodes           │   │  • genetic_profiles    │
│  • Agents         │   │  • Pathway nodes       │   │  • user_snps           │
│  • Trajectories   │   │  • Nutrient nodes      │   │  • lab_results         │
│  • Memory entries │   │  • Herb nodes          │   │  • symptoms            │
│  • Hyperbolic     │   │  • Condition nodes     │   │  • treatments          │
│                   │   │  • Publication nodes   │   │  • subscriptions       │
│  HNSW indices     │   │  • HAS_VARIANT edges   │   │  • orders              │
│  Vector search    │   │  • AFFECTS edges       │   │  • practitioners       │
│                   │   │  • TREATS edges        │   │                        │
│  Port: 5432       │   │  • INVOLVED_IN edges   │   │  HIPAA-compliant       │
│  Docker volume    │   │  • Hyperedges          │   │  Field encryption      │
│                   │   │                        │   │  Port: 5433            │
└───────────────────┘   └────────────────────────┘   └────────────────────────┘
        │                           │                           │
        │                           │                           │
        └───────────────────────────┴───────────────────────────┘
                                    │
                        ┌───────────┴───────────┐
                        │  Cross-DB Query Layer │
                        │  (User SNPs → Graph)  │
                        └───────────────────────┘
```

---

## Database 1: RuVector PostgreSQL (Operational)

**Purpose:** Claude Flow learning, vector embeddings, AI intelligence

**Docker Setup:** (Already configured in `./ruvector-postgres/docker-compose.yml`)

```yaml
services:
  postgres:
    image: ruvnet/ruvector-postgres:latest
    container_name: ruvector-postgres
    ports:
      - "5432:5432"
    volumes:
      - postgres_data:/var/lib/postgresql/data
```

**Tables:**
- `embeddings` - Vector storage with HNSW (150x-12,500x faster search)
- `patterns` - ReasoningBank learned patterns
- `agents` - Multi-agent memory coordination
- `trajectories` - SONA reinforcement learning
- `memory_entries` - Claude Flow memory
- `hyperbolic_embeddings` - Hierarchical data (Poincare ball)
- `graph_nodes` / `graph_edges` - GNN message passing

**Performance:** ~61µs latency, 16,400 QPS with HNSW

---

## Database 2: RuVector Graph (Knowledge Graph)

**Purpose:** Biomedical knowledge graph connecting genes, SNPs, pathways, nutrients, herbs, conditions

**Integration with existing RuVector stack:**

```typescript
import { CodeGraph } from 'ruvector/dist/core/graph-wrapper';

// Initialize knowledge graph with persistent storage
const knowledgeGraph = new CodeGraph({
  storagePath: './data/knowledge-graph.db',
  inMemory: false
});

// Create nodes with typed labels
const mthfrGene = knowledgeGraph.createNode('HGNC:7436', ['Gene'], {
  symbol: 'MTHFR',
  name: 'Methylenetetrahydrofolate reductase',
  chromosome: '1',
  ensembl_id: 'ENSG00000177000'
});

const rs1801133 = knowledgeGraph.createNode('rs1801133', ['SNP'], {
  rsid: 'rs1801133',
  gene_symbol: 'MTHFR',
  chromosome: '1',
  position: 11856378,
  clinical_significance: 'drug_response',
  common_name: 'MTHFR C677T'
});

// Create typed edges
knowledgeGraph.createEdge('HGNC:7436', 'rs1801133', 'HAS_VARIANT', {
  effect: 'reduces_activity',
  evidence_level: 'clinical'
});
```

### Node Types (Matching Original Neo4j Schema)

| Node Type | Count Target | Properties |
|-----------|--------------|------------|
| Gene | ~20,000 | symbol, name, chromosome, ensembl_id, uniprot_id |
| SNP | ~1.2M | rsid, chromosome, position, ref_allele, alt_allele, maf |
| Pathway | ~2,000 | name, description, reactome_id, kegg_id |
| Nutrient | ~1,000 | name, category, usda_id, recommended_daily |
| Herb | ~2,000 | name, latin_name, modality, tcmsp_id, properties |
| Condition | ~5,000 | name, icd10_code, mesh_id, category |
| Publication | ~500,000 | pmid, title, journal, year, evidence_level |
| Formula | ~500 | name, modality, category, indication |

### Edge Types

```typescript
// Gene relationships
(Gene)-[:HAS_VARIANT]->(SNP)
(Gene)-[:INVOLVED_IN {role: 'enzyme'}]->(Pathway)

// SNP relationships
(SNP)-[:AFFECTS {effect: 'reduces_activity', magnitude: 'moderate'}]->(Nutrient)
(SNP)-[:ASSOCIATED_WITH {odds_ratio: 1.5, p_value: 0.001}]->(Condition)

// Pathway relationships
(Pathway)-[:REQUIRES_COFACTOR]->(Nutrient)

// Herb relationships
(Herb)-[:TREATS {evidence_level: 'clinical'}]->(Condition)
(Herb)-[:INGREDIENT_OF {role: 'chief'}]->(Formula)

// Research relationships
(SNP|Herb|Nutrient)-[:HAS_RESEARCH]->(Publication)
```

### Hyperedges for Complex Relationships

RuVector Graph supports hyperedges - edges connecting multiple nodes simultaneously:

```typescript
// A hyperedge connecting multiple genes to a shared pathway
knowledgeGraph.createHyperedge(
  ['MTHFR', 'MTR', 'MTRR', 'SHMT1'],  // Multiple genes
  'METHYLATION_CYCLE_PARTICIPANTS',
  {
    pathway_id: 'REACT_R-HSA-196741',
    role: 'core_enzymes'
  }
);

// A hyperedge for TCM formula composition
knowledgeGraph.createHyperedge(
  ['huang_qi', 'dang_shen', 'bai_zhu', 'gan_cao'],  // Herbs in formula
  'FORMULA_COMPOSITION',
  {
    formula: 'Bu Zhong Yi Qi Tang',
    traditional_use: 'qi_tonifying'
  }
);
```

### Cypher Queries (Compatible Syntax)

```typescript
// Find nutrients affected by user's SNPs
const result = knowledgeGraph.cypher(`
  MATCH (snp:SNP)-[a:AFFECTS]->(nutrient:Nutrient)
  WHERE snp.rsid IN $user_rsids
  RETURN snp.rsid, snp.common_name, nutrient.name, a.effect
  ORDER BY a.magnitude DESC
`, { user_rsids: ['rs1801133', 'rs1801131'] });

// Find TCM herbs for a condition via pathways
const tcmResult = knowledgeGraph.cypher(`
  MATCH (condition:Condition {name: $condition})
        <-[:TREATS]-(herb:Herb {modality: 'tcm'})
        -[:HAS_RESEARCH]->(pub:Publication)
  WHERE pub.evidence_level IN ['clinical_trial', 'meta_analysis']
  RETURN herb.name, herb.properties, count(pub) as research_count
  ORDER BY research_count DESC
  LIMIT 10
`, { condition: 'ME/CFS' });
```

### Graph Algorithms

```typescript
// Find most connected genes (hub identification)
const pageRanks = knowledgeGraph.pageRank(20, 0.85);
const topGenes = [...pageRanks.entries()]
  .filter(([id]) => id.startsWith('HGNC:'))
  .sort((a, b) => b[1] - a[1])
  .slice(0, 10);

// Detect pathway communities
const communities = knowledgeGraph.communities();  // Louvain algorithm

// Find shortest path between gene and condition
const path = knowledgeGraph.shortestPath('MTHFR', 'COND_depression', 5);

// Get all paths (for showing multiple mechanisms)
const allPaths = knowledgeGraph.allPaths('MTHFR', 'COND_depression', 4, 10);
```

---

## Database 3: PostgreSQL (User Domain - HIPAA)

**Purpose:** User PII, genetic profiles, health data with HIPAA compliance

**Port:** 5433 (to avoid conflict with RuVector PostgreSQL on 5432)

**Tables:** (From existing schema in 45-DATA-MODEL.md)
- `users` - Account information
- `genetic_profiles` - Upload metadata, encrypted raw data
- `user_snps` - User's SNP data (links to knowledge graph via rsid)
- `lab_results` - Lab test results
- `symptoms` - User-reported symptoms
- `treatments` - User treatments/supplements
- `subscriptions` - Subscription status
- `orders` - Marketplace orders
- `practitioners` - Practitioner profiles

**Security:**
- Field-level encryption (AES-256) for genetic data
- Separate database server (isolation)
- Audit logging
- HIPAA BAA-ready

---

## Docker Configuration

### Updated docker-compose.yml

```yaml
# RuVector Unified Stack
# Three-database architecture: RuVector PostgreSQL + RuVector Graph + User PostgreSQL

services:
  # Database 1: RuVector PostgreSQL (Operational/Vector)
  ruvector-postgres:
    image: ruvnet/ruvector-postgres:latest
    container_name: ruvector-postgres
    environment:
      POSTGRES_USER: claude
      POSTGRES_PASSWORD: ${RUVECTOR_PG_PASSWORD:-claude-flow-dev}
      POSTGRES_DB: claude_flow
    ports:
      - "5432:5432"
    volumes:
      - ruvector_postgres_data:/var/lib/postgresql/data
      - ./scripts/init-ruvector.sql:/docker-entrypoint-initdb.d/01-init.sql
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U claude -d claude_flow"]
      interval: 5s
      timeout: 5s
      retries: 10
    command: >
      postgres
      -c work_mem=256MB
      -c maintenance_work_mem=512MB
    networks:
      - gene-network

  # Database 2: RuVector Graph (Knowledge Graph)
  # Note: RuVector Graph runs as embedded library, not separate container
  # Data stored in mounted volume, accessed via Node.js application
  # This service just provides a persistent volume mount point

  # Database 3: PostgreSQL (User Domain - HIPAA)
  user-postgres:
    image: postgres:15-alpine
    container_name: user-postgres
    environment:
      POSTGRES_USER: gene_app
      POSTGRES_PASSWORD: ${USER_PG_PASSWORD:-gene-app-dev}
      POSTGRES_DB: gene_users
    ports:
      - "5433:5432"
    volumes:
      - user_postgres_data:/var/lib/postgresql/data
      - ./scripts/init-users.sql:/docker-entrypoint-initdb.d/01-init.sql
    healthcheck:
      test: ["CMD-SHELL", "pg_isready -U gene_app -d gene_users"]
      interval: 5s
      timeout: 5s
      retries: 10
    networks:
      - gene-network

  # Optional: pgAdmin for database management
  pgadmin:
    image: dpage/pgadmin4:latest
    container_name: gene-pgadmin
    environment:
      PGADMIN_DEFAULT_EMAIL: admin@gene.local
      PGADMIN_DEFAULT_PASSWORD: admin
      PGADMIN_CONFIG_SERVER_MODE: 'False'
    ports:
      - "5050:80"
    depends_on:
      ruvector-postgres:
        condition: service_healthy
      user-postgres:
        condition: service_healthy
    profiles:
      - gui
    networks:
      - gene-network

volumes:
  ruvector_postgres_data:
    name: gene-ruvector-postgres
  user_postgres_data:
    name: gene-user-postgres
  knowledge_graph_data:
    name: gene-knowledge-graph

networks:
  gene-network:
    name: gene-network
```

---

## Cross-Database Query Patterns

### Pattern 1: User SNPs to Knowledge Graph Lookup

```typescript
// 1. Get user's SNPs from User PostgreSQL
const userSnps = await userDb.query(`
  SELECT rsid, genotype, risk_level
  FROM user_snps
  WHERE genetic_profile_id = $1
    AND risk_level IN ('heterozygous', 'homozygous_risk')
`, [profileId]);

// 2. Query knowledge graph for affected nutrients
const rsids = userSnps.map(s => s.rsid);
const affected = knowledgeGraph.cypher(`
  MATCH (snp:SNP)-[a:AFFECTS]->(nutrient:Nutrient)
  WHERE snp.rsid IN $rsids
  RETURN snp.rsid, nutrient.name, a.effect, a.magnitude
`, { rsids });

// 3. Combine results
const recommendations = userSnps.map(snp => ({
  ...snp,
  affectedNutrients: affected.rows
    .filter(r => r[0] === snp.rsid)
    .map(r => ({ name: r[1], effect: r[2], magnitude: r[3] }))
}));
```

### Pattern 2: Semantic Search with Graph Context

```typescript
// 1. User asks: "What supplements help with methylation if I have MTHFR?"
const queryEmbedding = await embedder.embed(userQuery);

// 2. Search RuVector PostgreSQL for similar patterns
const similarPatterns = await ruvectorDb.query(`
  SELECT content, metadata
  FROM embeddings
  ORDER BY vector <-> $1
  LIMIT 5
`, [queryEmbedding]);

// 3. Expand with graph relationships
const graphContext = knowledgeGraph.cypher(`
  MATCH (g:Gene {symbol: 'MTHFR'})-[:HAS_VARIANT]->(s:SNP)
        -[:AFFECTS]->(n:Nutrient)
  RETURN g.symbol, s.rsid, n.name, n.recommended_daily
`);

// 4. Combine for RAG context
const ragContext = {
  vectorResults: similarPatterns,
  graphContext: graphContext.rows,
  userSnps: userProfile.snps.filter(s => s.gene === 'MTHFR')
};
```

### Pattern 3: Pathway Visualization with User Overlay

```typescript
// 1. Get pathway structure from graph
const pathwayNodes = knowledgeGraph.cypher(`
  MATCH (p:Pathway {name: 'Methylation Cycle'})
        <-[:INVOLVED_IN]-(g:Gene)
        -[:HAS_VARIANT]->(s:SNP)
  RETURN g.symbol, s.rsid, s.clinical_significance
`);

// 2. Get user's genotypes for these SNPs
const rsids = pathwayNodes.rows.map(r => r[1]);
const userGenotypes = await userDb.query(`
  SELECT rsid, genotype, risk_level
  FROM user_snps
  WHERE genetic_profile_id = $1
    AND rsid = ANY($2)
`, [profileId, rsids]);

// 3. Merge for Cytoscape.js visualization
const cytoscapeData = buildCytoscapeGraph(pathwayNodes, userGenotypes);
```

---

## Data Migration Plan

### Phase 1: Knowledge Graph Population (Batch)

```typescript
import { CodeGraph } from 'ruvector/dist/core/graph-wrapper';
import { buildGraph } from 'ruvector/dist/core/graph-algorithms';

async function populateKnowledgeGraph() {
  const graph = new CodeGraph({
    storagePath: './data/knowledge-graph.db'
  });

  // 1. Import genes from HGNC
  console.log('Importing genes...');
  for await (const gene of hgncStream()) {
    graph.createNode(gene.hgnc_id, ['Gene'], gene);
  }

  // 2. Import SNPs from dbSNP/ClinVar
  console.log('Importing SNPs...');
  for await (const snp of dbSnpStream()) {
    graph.createNode(snp.rsid, ['SNP'], snp);
    if (snp.gene_id) {
      graph.createEdge(snp.gene_id, snp.rsid, 'HAS_VARIANT');
    }
  }

  // 3. Import pathways from Reactome/KEGG
  // 4. Import nutrients from USDA/FooDB
  // 5. Import herbs from TCMSP/IMPPAT
  // 6. Import conditions from MeSH/ICD-10
  // 7. Create relationship edges

  graph.save();
  console.log('Knowledge graph populated:', graph.stats());
}
```

### Phase 2: Incremental Updates

```typescript
// Weekly update job
async function updateKnowledgeGraph() {
  const graph = new CodeGraph({
    storagePath: './data/knowledge-graph.db'
  });
  graph.load();

  // Check for new publications
  const newPubs = await fetchNewPubMedArticles(lastUpdateDate);
  for (const pub of newPubs) {
    graph.createNode(`PMID:${pub.pmid}`, ['Publication'], pub);
    // Create HAS_RESEARCH edges based on MeSH terms
  }

  graph.save();
}
```

---

## Performance Comparison

| Operation | RuVector Graph | Neo4j | Improvement |
|-----------|---------------|-------|-------------|
| Node insert (batch) | 131K ops/sec | ~50K ops/sec | 2.6x faster |
| Edge traversal | <0.1ms | ~1ms | 10x faster |
| Cypher query | Compatible | Native | Equivalent |
| Vector+Graph hybrid | Native | Plugin required | Simpler |
| Memory per node | ~50 bytes | ~100 bytes | 2x efficient |

---

## Implementation Steps

### Phase 1: Infrastructure Setup (Week 1)

1. [x] Update docker-compose.yml with three-database architecture
2. [ ] Create init-users.sql for User PostgreSQL schema
3. [ ] Initialize RuVector Graph storage directory
4. [ ] Configure environment variables

### Phase 2: Knowledge Graph Population (Weeks 2-3)

1. [ ] Create data import scripts for each source
2. [ ] Import genes from HGNC (~20K nodes)
3. [ ] Import SNPs from dbSNP/ClinVar (~1.2M nodes)
4. [ ] Import pathways from Reactome/KEGG (~2K nodes)
5. [ ] Import nutrients, herbs, conditions
6. [ ] Create all relationship edges
7. [ ] Validate with sample queries

### Phase 3: Application Integration (Week 4)

1. [ ] Create database connection utilities
2. [ ] Implement cross-database query patterns
3. [ ] Build Cytoscape.js visualization data layer
4. [ ] Integrate with RAG pipeline

### Phase 4: Testing & Optimization (Week 5)

1. [ ] Performance benchmarks
2. [ ] Query optimization
3. [ ] Load testing with realistic data
4. [ ] Documentation

---

## Unified RuVector Stack Benefits

1. **Single Vendor Stack** - All vector and graph operations through RuVector
2. **Native Performance** - Rust bindings, 10x faster than WASM alternatives
3. **Vector-Graph Hybrid** - Seamless integration of embeddings with graph traversal
4. **Cypher Compatibility** - Same query syntax as Neo4j for easy migration
5. **Hypergraph Support** - Model complex n-ary relationships (TCM formulas, pathways)
6. **Reduced Complexity** - No Neo4j infrastructure to manage
7. **Cost Savings** - No Neo4j licensing, simpler deployment

---

## Hive Discussion: Graph Database Comparison (January 2026)

### Verified FREE Graph Databases (Apache/MIT/MPL Only)

| Database | License | Query Lang | Hypergraph | Vector | Rust Client | Best For |
|----------|---------|------------|------------|--------|-------------|----------|
| **RuVector Graph** | Apache 2.0 | Cypher | **YES** | **HNSW** | Native | Unified vector+graph |
| **Apache AGE** | Apache 2.0 | Cypher+SQL | No | No | `apache_age` | PostgreSQL shops |
| **NebulaGraph** | Apache 2.0 | nGQL/GQL | No | No | `nebula-rust` | Massive scale |
| **TypeDB** | MPL 2.0 | TypeQL | **YES** | No | `typedb-driver` | Knowledge graphs |
| **CozoDB** | MPL 2.0 | Datalog | No | **HNSW** | Native | Graph+Vector |
| **Oxigraph** | MIT/Apache | SPARQL | No | No | Native | RDF stores |
| **IndraDB** | MPL 2.0 | gRPC API | No | No | Native | TAO-style |
| **JanusGraph** | Apache 2.0 | Gremlin | No | No | `gremlin-client` | Distributed |
| **agdb** | Apache 2.0 | Builder | No | No | Native | Embedded ACID |

### Eliminated (NOT Free/Commercial)

- ~~Neo4j Enterprise~~ - $65K+/year commercial
- ~~TigerGraph~~ - Commercial, expensive
- ~~ArangoDB~~ - BSL license, 100GB commercial limit
- ~~Memgraph~~ - BSL license
- ~~SurrealDB~~ - BSL license
- ~~FalkorDB~~ - SSPLv1 (restrictive)
- ~~Amazon Neptune~~ - AWS commercial

---

### RuVector Deep-Dive Findings

#### Storage Backends (6 Options - Updated v3.0)

| Backend | Type | Crate | Use Case | Best For |
|---------|------|-------|----------|----------|
| **redb ^2.1** | Embedded KV | Primary | Persistent storage | Single-node, CLI, IoT |
| **memmap2 ^0.9** | Memory-mapped | Large datasets | Vectors > RAM | Millions of vectors |
| **PostgreSQL** | Server | ruvector-postgres | 77+ SQL functions | Enterprise, SQL teams |
| **In-memory** | Transient | Hot data | Caching | Testing, Lambda |
| **WASM/IndexedDB** | Browser | @ruvector/wasm | Offline-first | Privacy-first apps |
| **Distributed (Raft)** | Cluster | ruvector_raft | High availability | Mission-critical |

#### Storage Backend Selection Guide

| Requirement | Recommended Backend |
|-------------|---------------------|
| Single-node embedded | redb + memmap2 |
| Browser/client-side | WASM + IndexedDB |
| Enterprise/SQL needed | PostgreSQL (ruvector-postgres) |
| High availability | Distributed (Raft) |
| Maximum performance | In-memory |
| Large datasets (>RAM) | memmap2 |
| Testing/development | In-memory |

#### Tiered Storage Architecture

RuVector uses automatic tiered storage requiring **no configuration**:

| Tier | Access Pattern | Storage Type | Precision |
|------|----------------|--------------|-----------|
| **Hot** | Frequently accessed | In-memory cache | Full precision |
| **Warm** | Moderate access | Memory-mapped files | Efficient |
| **Cold** | Rarely accessed | Compressed on disk | Background optimized |

#### Confirmed Rust Dependencies

- **petgraph ^0.6** - Graph data structures
- **hnsw_rs ^0.3** - HNSW approximate nearest neighbor
- **simsimd ^5.9** - SIMD distance calculations
- **rayon ^1.10** - Data parallelism
- **redb ^2.1** - Embedded storage

#### Cypher Support (Partial)

| Supported | Not Supported |
|-----------|---------------|
| MATCH, OPTIONAL MATCH | Subqueries |
| CREATE, MERGE, SET | FOREACH loops |
| DELETE, DETACH DELETE | UNWIND |
| RETURN, WITH, WHERE | BEGIN/COMMIT |
| ORDER BY, SKIP, LIMIT | APOC procedures |
| Variable-length paths | Neo4j 5.x advanced |

#### Graph Algorithm Gaps (CRITICAL)

| Algorithm | RuVector | Required For | Alternative |
|-----------|----------|--------------|-------------|
| **PageRank** | **NO** | Hub gene identification | Use petgraph directly |
| **Shortest Path** | **NO** | Pathway traversal | Use petgraph directly |
| **Betweenness Centrality** | **NO** | Key connector nodes | Use petgraph directly |
| **Community Detection** | **NO** | Pathway clustering | Use petgraph directly |
| **Connected Components** | **NO** | Network analysis | Use petgraph directly |

**What RuVector HAS:**
- HNSW vector search (16,400 QPS)
- Min-Cut algorithms (Karger-Stein, Stoer-Wagner)
- GNN layers (GCN, GraphSAGE, GAT)
- 39 attention mechanisms
- Hyperedge support (N-ary relationships)

---

### Scalability Analysis

#### Memory Requirements for 1.7M Nodes

| Configuration | Memory | Notes |
|---------------|--------|-------|
| Vectors (PQ8 compressed) | ~340 MB | 8x compression |
| Vectors (f32 uncompressed) | ~2.6 GB | Full precision |
| HNSW index overhead | ~5 GB | 2x vector size |
| Graph structure | ~500 MB | Edges, properties |
| **Production Total** | **10-12 GB** | Recommended |

#### Server Sizing

| Tier | Nodes | RAM | CPU | Use Case |
|------|-------|-----|-----|----------|
| Development | <100K | 8 GB | 4 cores | Local dev |
| Staging | ~1.7M | 16 GB | 8 cores | Testing |
| **Production** | ~1.7M | **32 GB** | **16 cores** | Live traffic |
| Growth (2yr) | ~5M | 64 GB | 32 cores | Scale target |

#### Scaling Ceiling

| Scale | Feasibility | Architecture |
|-------|-------------|--------------|
| 1.7M nodes | **YES** | Single embedded instance |
| 5M nodes | **YES** | Larger VM, 64GB RAM |
| 10M nodes | **CAUTION** | Consider sharding |
| 100M+ nodes | **NO** | Requires distributed |

---

### QLever Analysis: Biomedical SPARQL Engine (New in v3.0)

#### Overview

**QLever** is a high-performance SPARQL engine (Apache 2.0) developed at University of Freiburg. It scales to **1 trillion+ triples** on commodity hardware.

| Property | Value |
|----------|-------|
| License | Apache 2.0 (fully open source) |
| Max Scale | 1T+ triples |
| RAM (7B triples) | ~40 GB |
| Index Speed | ~1hr/1B triples |
| Vector Search | **NO** - Critical gap |
| Docker | Yes (adfreiburg/qlever) |

#### Proven Biomedical Deployments

- **UniProt** - Protein sequence database (billions of triples)
- **PubChem** - Chemical compound database
- **dblp** - Official SPARQL endpoint

#### Verdict for Gene Platform

**Partial Fit** - QLever excels at SPARQL queries on biomedical RDF (gene ontologies, pathway relationships) but **cannot provide vector/embedding search**.

**Recommendation**: Use QLever in **hybrid architecture** alongside RuVector:

```
┌─────────────────┐          ┌─────────────────────────┐
│    QLever       │          │      RuVector           │
│  (SPARQL/RDF)   │          │   (Vector/HNSW)         │
├─────────────────┤          ├─────────────────────────┤
│ Gene ontologies │          │ Gene embeddings         │
│ Pathway triples │          │ SNP similarity vectors  │
│ Federated data  │          │ Semantic search         │
└─────────────────┘          └─────────────────────────┘
```

---

### Critical Analysis: Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Missing graph algorithms | HIGH | HIGH | Integrate petgraph directly |
| 62% documentation coverage | HIGH | MEDIUM | Budget reverse-engineering time |
| Version 0.1.x stability | MEDIUM | HIGH | Extensive testing, escape plan |
| No APOC procedures | HIGH | MEDIUM | Build custom functions |
| Vendor dependency | LOW | MEDIUM | Cypher enables migration path |
| **Issue #103: SIMD garbled output** | HIGH | CRITICAL | Verify fixed before adoption |
| **Issue #110: No ARM64 support** | HIGH | MEDIUM | Use x86_64 or alternative |
| **Single maintainer (ruvnet)** | MEDIUM | HIGH | Fork repository, build expertise |

#### Proof-of-Concept Requirements (Before Commit)

1. **2-week POC** with real biomedical data (10M edges minimum)
2. **Benchmark** 6-hop traversals, PageRank (via petgraph), vector+graph hybrid
3. **Verify** backup/restore, memory under load, query profiling
4. **Document** escape plan to Apache AGE or CozoDB

---

### Decision Framework (Updated by Hive-Mind Analysis)

#### Weighted Criteria (Revised v3.0)

**Scoring Methodology**: 1 (Poor) to 5 (Excellent)

| Criterion | Weight | RuVector | CozoDB | Apache AGE | TypeDB |
|-----------|--------|----------|--------|------------|--------|
| **Maturity** | 25% | 1 | 2 | 4 | 5 |
| **Vector Support** | 20% | 5 | 4 | 3* | 2 |
| **Graph Queries** | 15% | 4 | 4 | 4 | 5 |
| **Community** | 15% | 1 | 3 | 5 | 4 |
| **Documentation** | 10% | 2 | 3 | 4 | 4 |
| **Operational** | 10% | 1 | 2 | 5 | 4 |
| **Innovation** | 5% | 5 | 4 | 2 | 3 |
| **WEIGHTED SCORE** | 100% | **2.30** | **3.05** | **3.95** | **3.95** |

*Apache AGE requires pgvector separately for vector support

#### Critical Analysis Findings

| Risk Factor | RuVector Assessment | Severity |
|-------------|---------------------|----------|
| Version 0.1.x (pre-1.0) | No stability guarantees | HIGH |
| Single primary maintainer | Bus factor critical | HIGH |
| Issue #103: SIMD garbled output | Core inference broken | HIGH |
| Issue #110: No ARM64 binaries | No Apple Silicon | MEDIUM |
| 62% documentation coverage | Reverse-engineering needed | MEDIUM |

#### Revised Recommendation (v3.0)

**AVOID RuVector for production** (90% confidence) - Too immature, single maintainer, critical bugs

**Primary Recommendations:**
1. **TypeDB** (3.95) - Best for enterprise, production-critical, corporate support
2. **Apache AGE** (3.95) - Best for PostgreSQL shops, operational simplicity

**Secondary Options:**
3. **CozoDB** (3.05) - Best for embedded, experimental (Datalog + HNSW, MPL 2.0)
4. **RuVector** (2.30) - Only for cutting-edge AI features (accept HIGH risk)

**If RuVector is chosen anyway (Mandatory Mitigations):**
1. Fork repository immediately
2. Dedicate engineer to understand codebase
3. Minimum 3-week POC before production commitment
4. Document escape plan to Apache AGE or CozoDB
5. Implement custom health checks
6. Daily exports to portable format

---

### Implementation Plan Update

#### Phase 1: Proof-of-Concept (2 weeks)
- [ ] Deploy RuVector Graph with 1M test nodes
- [ ] Benchmark 6-hop pathway traversals
- [ ] Implement PageRank via petgraph wrapper
- [ ] Test hyperedge queries for TCM formulas
- [ ] Validate backup/restore procedures

#### Phase 2: Algorithm Integration (1 week)
- [ ] Create petgraph wrapper for missing algorithms
- [ ] Implement PageRank, shortest path, community detection
- [ ] Expose via Cypher UDFs or API layer

#### Phase 3: Production Setup (1 week)
- [ ] Configure 32GB production instance
- [ ] Setup monitoring and alerting
- [ ] Document escape plan procedures

---

## Escape Plan: If RuVector Fails (New in v3.0)

### Migration Path Options (Ranked by Effort) - FREE OPTIONS ONLY

| Path | Difficulty | Time Estimate | Query Rewrite | License |
|------|------------|---------------|---------------|---------|
| RuVector → Apache AGE | EASY | 2-4 weeks | Minor (OpenCypher) | Apache 2.0 |
| RuVector → CozoDB | HARD | 4-8 weeks | Complete (Datalog) | MPL 2.0 |
| RuVector → TypeDB | HARD | 6-10 weeks | Complete (TypeQL) | MPL 2.0 |

**NOT viable (commercial/restrictive):**
- ~~Neo4j~~ - $65K+/year Enterprise, Community has restrictions
- ~~TigerGraph~~ - Commercial
- ~~ArangoDB~~ - BSL license, 100GB commercial limit
- ~~Memgraph~~ - BSL license
- ~~SurrealDB~~ - BSL license
- ~~FalkorDB~~ - SSPLv1 (restrictive)

### Recommended Primary Escape: Apache AGE + pgvector

```cypher
-- RuVector (Cypher)
MATCH (p:Person)-[:KNOWS]->(f:Person)
WHERE p.name = 'Alice'
RETURN f.name

-- Apache AGE (OpenCypher) - Nearly identical wrapper
SELECT * FROM cypher('graph_name', $$
MATCH (p:Person)-[:KNOWS]->(f:Person)
WHERE p.name = 'Alice'
RETURN f.name
$$) AS (friend_name agtype);
```

**Migration Effort**: 80% of queries need minimal modification.

### Pre-Migration Steps (Document Now)

1. Document all schema relationships
2. Export vectors to standard format (JSONL, Parquet)
3. Maintain query log for migration testing
4. Create test fixtures from production data

---

## Open Questions (Updated v3.0)

- [x] ~~Is RuVector PostgreSQL sufficient for graph operations?~~ **No - use RuVector Graph for knowledge graph**
- [x] ~~How to integrate petgraph for missing algorithms?~~ **Use petgraph directly; RuVector lacks PageRank, centrality**
- [ ] **NEW (v3.0):** Given 2.30 score, should we adopt CozoDB or Apache AGE instead?
- [ ] **NEW (v3.0):** Is QLever hybrid architecture viable for biomedical SPARQL queries?
- [x] ~~How large will the dataset grow?~~ **~1.7M nodes, within capacity with 32GB RAM**
- [x] ~~What query patterns are critical?~~ **Cypher for relationships, hybrid vector+graph for RAG**
- [x] ~~Should we consider CozoDB as alternative?~~ **Yes - scored 3.05, better than RuVector's 2.30**

---

## Files to Create/Update

- `./docker-compose.yml` - Updated three-database configuration
- `./scripts/init-users.sql` - User PostgreSQL schema
- `./src/lib/knowledge-graph.ts` - RuVector Graph wrapper
- `./src/lib/cross-db-queries.ts` - Cross-database query utilities
- `./scripts/populate-knowledge-graph.ts` - Data import scripts

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | Jan 2026 | Claude | Initial plan based on research |
| 2.0 | Jan 2026 | Claude | Replace Neo4j with RuVector Graph |
| 3.0 | Jan 2026 | Claude | Hive-mind analysis (7 agents): Critical findings on RuVector maturity (2.30 score), expanded storage backends (6 options), added QLever SPARQL analysis, escape plan, revised recommendations favoring TypeDB/Apache AGE |

---

## Hive-Mind Research Documentation (v3.0)

**Research Location:** `docs/research/graphs/`

| Document | Lines | Content |
|----------|-------|---------|
| 01-RUVECTOR-DEEP-DIVE.md | 731 | Cargo.toml analysis, Cypher support, algorithms |
| 02-FREE-GRAPH-DATABASES.md | 560 | 8 free alternatives compared |
| 03-RUST-GRAPH-IMPLEMENTATIONS.md | 1,015 | petgraph, CozoDB, agdb analysis |
| 04-SCALABILITY-ANALYSIS.md | 904 | Memory calculations, server sizing |
| 05-CRITICAL-ANALYSIS.md | 480 | Risk assessment, weighted scoring |
| 06-RUVECTOR-STORAGE-BACKENDS.md | 630 | 6 storage backends documented |
| 07-QLEVER-ANALYSIS.md | 473 | SPARQL engine for biomedical data |
| **Total** | **4,793** | **172 KB persisted to disk** |
