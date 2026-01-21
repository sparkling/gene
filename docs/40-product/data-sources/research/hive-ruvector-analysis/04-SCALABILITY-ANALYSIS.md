# Scalability & Performance Analysis for RuVector-Based Knowledge Graph

**Analysis Date:** 2026-01-21
**Analyst:** System Architecture Designer (Hive-Mind Agent)
**Target System:** Gene-Disease Knowledge Graph with Vector Search

---

## Target Scale from Data Model

Based on `/docs/40-product/45-DATA-MODEL.md`:

| Entity Type | Count | Notes |
|-------------|-------|-------|
| Genes | ~20,000 | Human protein-coding genes |
| SNPs | ~1,200,000 | Single nucleotide polymorphisms |
| Publications | ~500,000 | PubMed references |
| Pathways | ~2,000 | KEGG, Reactome pathways |
| Conditions | ~5,000 | Diseases, phenotypes |
| **Total Nodes** | **~1,727,000** | |
| **Estimated Edges** | **~10,000,000+** | Gene-SNP, Gene-Pathway, etc. |

---

## 1. Memory Requirements Calculation

### Base Data Size Calculations

**Vector Embedding Storage:**
- Embedding dimension: 384 (MiniLM-L6-v2 standard)
- Bytes per float32: 4 bytes
- Vector size per node: 384 x 4 = **1,536 bytes (1.5 KB)**

**Node Counts & Embedding Requirements:**

| Entity | Count | Vector Size | Total Vectors |
|--------|-------|-------------|---------------|
| Genes | 20,000 | 1.5 KB | 30 MB |
| SNPs | 1,200,000 | 1.5 KB | 1,800 MB |
| Publications | 500,000 | 1.5 KB | 750 MB |
| Pathways | 2,000 | 1.5 KB | 3 MB |
| Conditions | 5,000 | 1.5 KB | 7.5 MB |
| **Total Vectors** | **1,727,000** | | **2,590 MB (~2.5 GB)** |

---

### HNSW Index Overhead

HNSW (Hierarchical Navigable Small World) indices typically require:
- **M parameter** (connections per node): 16-32 typical
- **ef_construction**: 200 typical
- **Overhead multiplier**: 1.5x - 2x vector storage

**HNSW Index Size Calculation:**
```
HNSW overhead = vector_size × nodes × M × 8 bytes (link storage)
             = 384 × 1,727,000 × 16 × 8
             = ~84 GB for pure HNSW graph structure

Practical estimate with compression: ~3-5 GB additional
```

**Total Vector + Index Storage:**
- Vectors: 2.5 GB
- HNSW Index: 4 GB (practical estimate with optimizations)
- **Subtotal: ~6.5 GB**

---

### Graph Structure Overhead

**Edge Storage (10M+ edges):**
- Edge ID: 8 bytes
- Source node ID: 8 bytes
- Target node ID: 8 bytes
- Edge type: 4 bytes (enum)
- Properties (avg): ~100 bytes
- **Per edge: ~128 bytes**

**Edge Storage Total:**
```
10,000,000 edges × 128 bytes = 1.28 GB
```

**Node Properties (beyond embeddings):**
- Average node properties: ~500 bytes
- 1,727,000 nodes × 500 bytes = ~864 MB

**Graph Indexes:**
- Node ID index: ~50 MB
- Label index: ~100 MB
- Property indexes (symbol, rsid, name): ~500 MB

**Subtotal Graph: ~2.8 GB**

---

### Total Memory Requirements by Database Type

#### Option 1: RuVector (Embedded TypeScript)

| Component | Size | Notes |
|-----------|------|-------|
| Vector embeddings | 2.5 GB | 384-dim float32 |
| HNSW index | 4.0 GB | M=16, optimized |
| Graph structure | 2.8 GB | Edges + properties |
| Working memory | 2.0 GB | Query buffers, caches |
| **Total** | **~11.3 GB** | Single-node requirement |

#### Option 2: Kuzu (Embedded C++)

| Component | Size | Notes |
|-----------|------|-------|
| Vector embeddings | 2.5 GB | Native arrays |
| Graph structure | 1.8 GB | Columnar compression |
| Property storage | 1.5 GB | Compressed |
| Working memory | 1.5 GB | Buffer pools |
| **Total** | **~7.3 GB** | More efficient storage |

*Note: Kuzu doesn't have native HNSW - would need external vector search*

#### Option 3: Neo4j + pgvector (Separate Systems)

| Component | Size | Notes |
|-----------|------|-------|
| Neo4j heap | 4.0 GB | Minimum for this scale |
| Neo4j page cache | 4.0 GB | For graph traversals |
| pgvector storage | 2.5 GB | Vectors |
| pgvector HNSW | 4.0 GB | Index |
| PostgreSQL shared_buffers | 2.0 GB | |
| **Total** | **~16.5 GB** | Split across 2 processes |

#### Option 4: SurrealDB

| Component | Size | Notes |
|-----------|------|-------|
| Document storage | 3.0 GB | JSON overhead |
| Graph edges | 2.0 GB | |
| Indexes | 2.5 GB | Multiple index types |
| Vector store | 2.5 GB | |
| HNSW index | 4.0 GB | |
| **Total** | **~14 GB** | All-in-one |

---

### Memory Summary Table

| Database | Minimum RAM | Recommended | 2-Year Growth |
|----------|-------------|-------------|---------------|
| RuVector | 16 GB | 24 GB | 32 GB |
| Kuzu | 12 GB | 16 GB | 24 GB |
| Neo4j + pgvector | 24 GB | 32 GB | 48 GB |
| SurrealDB | 20 GB | 28 GB | 40 GB |

---

## 2. Server Sizing Recommendations

### Development Environment

**Purpose:** Local development, testing, prototyping

| Database | CPU | RAM | Storage | Cost/Month |
|----------|-----|-----|---------|------------|
| RuVector | 2 vCPU | 8 GB | 50 GB SSD | N/A (local) |
| Kuzu | 2 vCPU | 4 GB | 30 GB SSD | N/A (local) |
| Neo4j + pgvector | 4 vCPU | 16 GB | 100 GB SSD | N/A (local) |
| SurrealDB | 2 vCPU | 8 GB | 50 GB SSD | N/A (local) |

**Dev Notes:**
- Use subset of data (~10% = 170K nodes) for local dev
- Docker Compose recommended for Neo4j + pgvector stack
- Embedded options (RuVector, Kuzu) simplest to set up

**Recommended Dev Machine Specs:**
- MacBook Pro M2/M3 with 16GB RAM minimum
- Or Linux workstation with 16GB+ RAM
- VS Code or similar IDE

---

### Staging Environment

**Purpose:** Integration testing, performance testing, demo

| Database | Instance Type (AWS) | vCPU | RAM | Storage | Est. Cost/Month |
|----------|---------------------|------|-----|---------|-----------------|
| RuVector | r6g.xlarge | 4 | 32 GB | 200 GB gp3 | ~$180 |
| Kuzu | r6g.large | 2 | 16 GB | 100 GB gp3 | ~$100 |
| Neo4j | r6g.2xlarge | 8 | 64 GB | 300 GB gp3 | ~$380 |
| pgvector | r6g.xlarge | 4 | 32 GB | 200 GB gp3 | ~$180 |
| SurrealDB | r6g.xlarge | 4 | 32 GB | 200 GB gp3 | ~$180 |

**Staging Stack Cost Comparison:**

| Option | Monthly Cost | Notes |
|--------|--------------|-------|
| RuVector (single) | ~$180 | Simplest |
| Kuzu + usearch | ~$150 | Most economical |
| Neo4j + pgvector | ~$560 | Two instances required |
| SurrealDB | ~$180 | All-in-one |

---

### Production Environment

**Purpose:** Live system with full dataset, high availability considerations

#### Option A: RuVector Production

| Component | Instance | vCPU | RAM | Storage | HA Config |
|-----------|----------|------|-----|---------|-----------|
| Primary | r6g.2xlarge | 8 | 64 GB | 500 GB gp3 | Active |
| Read Replica | r6g.xlarge | 4 | 32 GB | 500 GB gp3 | Passive (optional) |
| **Total** | | 12 | 96 GB | 1 TB | |

**Monthly Cost:** ~$500-750

#### Option B: Kuzu + External Vector Search

| Component | Instance | vCPU | RAM | Storage | Purpose |
|-----------|----------|------|-----|---------|---------|
| Kuzu | r6g.xlarge | 4 | 32 GB | 300 GB | Graph |
| Qdrant/Weaviate | r6g.xlarge | 4 | 32 GB | 300 GB | Vectors |
| **Total** | | 8 | 64 GB | 600 GB | |

**Monthly Cost:** ~$400

#### Option C: Neo4j + pgvector Production

| Component | Instance | vCPU | RAM | Storage | Config |
|-----------|----------|------|-----|---------|--------|
| Neo4j | r6g.4xlarge | 16 | 128 GB | 500 GB | Primary |
| Neo4j Read | r6g.2xlarge | 8 | 64 GB | 500 GB | Replica |
| pgvector (RDS) | db.r6g.2xlarge | 8 | 64 GB | 500 GB | Multi-AZ |
| **Total** | | 32 | 256 GB | 1.5 TB | |

**Monthly Cost:** ~$2,000-2,500

#### Option D: SurrealDB Production

| Component | Instance | vCPU | RAM | Storage | Config |
|-----------|----------|------|-----|---------|--------|
| Primary | r6g.2xlarge | 8 | 64 GB | 500 GB | Active |
| Secondary | r6g.xlarge | 4 | 32 GB | 500 GB | Standby |
| **Total** | | 12 | 96 GB | 1 TB | |

**Monthly Cost:** ~$500-700

---

### Production Cost Summary

| Option | Monthly Infra | Operational Complexity | Vendor Lock-in |
|--------|---------------|------------------------|----------------|
| RuVector | ~$500-750 | Low | Low |
| Kuzu + Qdrant | ~$400 | Medium | Low |
| Neo4j + pgvector | ~$2,000-2,500 | High | High (Neo4j) |
| SurrealDB | ~$500-700 | Low | Medium |

---

### 2-Year Growth Projections

**Assumptions:**
- Data growth: 50% annually
- User base growth: 100% annually
- Query volume: 200% annually

**Year 1 -> Year 2 Data Scale:**

| Metric | Year 1 | Year 2 | Growth |
|--------|--------|--------|--------|
| Nodes | 1.7M | 3.8M | +124% |
| Edges | 10M | 25M | +150% |
| Vectors | 1.7M | 3.8M | +124% |
| Daily queries | 100K | 500K | +400% |

**Year 2 Infrastructure Requirements:**

| Option | RAM | Storage | Monthly Cost |
|--------|-----|---------|--------------|
| RuVector | 64-96 GB | 1.5 TB | ~$1,000-1,500 |
| Kuzu + Vector DB | 48-64 GB | 1.2 TB | ~$800 |
| Neo4j + pgvector | 192-256 GB | 3 TB | ~$4,000-5,000 |
| SurrealDB | 64-96 GB | 1.5 TB | ~$1,000-1,500 |

---

## 3. Query Performance Expectations

### Critical Query Patterns for Gene-Disease Knowledge Graph

Based on the application requirements, these are the key query patterns:

---

### Query Type 1: 6-Hop Pathway Traversals

**Use Case:** Find all connections between a gene and distant conditions through multiple pathways

```cypher
MATCH path = (gene:Gene {symbol: 'MTHFR'})
    -[:HAS_VARIANT]->(snp:SNP)
    -[:ASSOCIATED_WITH]->(condition1:Condition)
    -[:TREATED_BY]->(herb:Herb)
    -[:AFFECTS]->(pathway:Pathway)
    -[:INVOLVES]->(gene2:Gene)
WHERE gene2 <> gene
RETURN path
LIMIT 100
```

**Expected Performance by Database:**

| Database | Cold Query | Warm Query | Notes |
|----------|------------|------------|-------|
| RuVector | 500-2000ms | 100-500ms | In-memory graph |
| Kuzu | 200-800ms | 50-200ms | Columnar optimization |
| Neo4j | 100-500ms | 20-100ms | Mature query optimizer |
| SurrealDB | 800-3000ms | 200-800ms | Document overhead |

**Analysis:**
- Neo4j excels here due to decades of graph optimization
- Kuzu competitive due to columnar storage for scans
- RuVector adequate but not optimized for deep traversals
- SurrealDB struggles with complex multi-hop queries

---

### Query Type 2: Vector Similarity + Graph Filter

**Use Case:** Find similar research papers that mention specific genes

```sql
-- Step 1: Vector search
SELECT id, embedding <=> query_embedding AS distance
FROM publications
WHERE distance < 0.3
ORDER BY distance
LIMIT 1000;

-- Step 2: Graph filter
MATCH (pub:Publication)-[:MENTIONS]->(gene:Gene)
WHERE pub.id IN $vector_results
  AND gene.symbol IN ['MTHFR', 'MTR', 'MTRR']
RETURN pub, gene
```

**Expected Performance by Database:**

| Database | Vector Search | Graph Filter | Combined |
|----------|---------------|--------------|----------|
| RuVector | 5-20ms (HNSW) | 50-200ms | 55-220ms |
| Kuzu + Qdrant | 2-10ms | 20-80ms | 30-100ms (2 RTT) |
| pgvector + Neo4j | 10-50ms | 20-100ms | 50-180ms (2 RTT) |
| SurrealDB | 20-100ms | 100-400ms | 120-500ms |

**Analysis:**
- Native integration (RuVector, SurrealDB) avoids network round-trips
- Kuzu + dedicated vector DB fastest but requires orchestration
- pgvector + Neo4j mature but operationally complex
- RuVector good balance of speed and simplicity

---

### Query Type 3: PageRank on 1.7M Nodes

**Use Case:** Identify most influential genes/conditions in the knowledge graph

```javascript
// Find hub genes
const pageRanks = await graph.pageRank({
  maxIterations: 20,
  dampingFactor: 0.85,
  nodeLabels: ['Gene'],
  relationshipTypes: ['HAS_VARIANT', 'INVOLVED_IN', 'AFFECTS']
});
```

**Expected Performance by Database:**

| Database | Time (1.7M nodes) | Memory Peak | Notes |
|----------|-------------------|-------------|-------|
| RuVector | 30-60 seconds | +2 GB | TypeScript implementation |
| Kuzu | 5-15 seconds | +1 GB | Native C++ algorithms |
| Neo4j GDS | 2-10 seconds | +4 GB | Highly optimized library |
| SurrealDB | 60-180 seconds | +3 GB | Not optimized for analytics |

**Analysis:**
- Neo4j Graph Data Science library is industry-leading for graph algorithms
- Kuzu's C++ core provides excellent algorithm performance
- RuVector functional but slower due to TypeScript
- SurrealDB not designed for heavy graph analytics

---

### Query Type 4: Batch Import Performance

**Use Case:** Initial data load and periodic updates

**Test Scenario:** Import 100,000 nodes with 384-dim embeddings + 500,000 edges

| Database | Node Import | Edge Import | Total Time | Notes |
|----------|-------------|-------------|------------|-------|
| RuVector | 10-20 min | 5-10 min | 15-30 min | Batch APIs available |
| Kuzu | 2-5 min | 1-3 min | 3-8 min | CSV bulk loader |
| Neo4j | 5-15 min | 3-8 min | 8-23 min | neo4j-admin import |
| SurrealDB | 15-30 min | 10-20 min | 25-50 min | Transaction overhead |

**Full Dataset Load (1.7M nodes, 10M edges):**

| Database | Estimated Time | Notes |
|----------|----------------|-------|
| RuVector | 4-8 hours | Streaming recommended |
| Kuzu | 1-2 hours | Fastest bulk loader |
| Neo4j | 2-4 hours | Use admin import tool |
| SurrealDB | 8-16 hours | Consider chunked loads |

---

### Query Performance Summary Matrix

| Query Type | RuVector | Kuzu | Neo4j | SurrealDB | Winner |
|------------|----------|------|-------|-----------|--------|
| 6-hop traversal | Good | Better | Best | Fair | Neo4j |
| Vector + Graph | Better | Best* | Good | Fair | Kuzu* |
| PageRank | Fair | Better | Best | Poor | Neo4j |
| Batch import | Good | Best | Better | Fair | Kuzu |
| Simple lookups | Best | Better | Good | Better | RuVector |
| Concurrent writes | Good | Fair | Better | Best | SurrealDB |

*Requires external vector DB

---

### Performance Optimization Strategies

#### For RuVector:
1. Enable HNSW indexing with M=16, ef=200
2. Use batch operations for writes
3. Implement query caching for repeated patterns
4. Consider read replicas for query scaling

#### For Kuzu:
1. Design schema for columnar access patterns
2. Use bulk loader for initial import
3. Integrate usearch/Qdrant for vector operations
4. Leverage parallel query execution

#### For Neo4j:
1. Tune heap and page cache sizes
2. Use GDS library for algorithms
3. Create appropriate indexes on frequently queried properties
4. Consider Enterprise for read replicas

#### For SurrealDB:
1. Use live queries for real-time updates
2. Implement proper indexing strategy
3. Consider document design for common access patterns
4. Use transactions judiciously

---

## 4. Scaling Ceiling Analysis

### Single-Node Limits

Understanding when each database option hits its ceiling is critical for long-term planning.

---

### RuVector Scaling Limits

**Architecture:** Embedded TypeScript/Node.js with better-sqlite3 backend

| Constraint | Limit | Impact |
|------------|-------|--------|
| Node.js heap | ~4 GB default, 8 GB practical | All indices must fit |
| SQLite file size | 281 TB theoretical, ~1 TB practical | Disk I/O bottleneck |
| Single-threaded writes | 1 thread | Write throughput capped |
| HNSW memory | ~2x vector size | Memory dominant factor |

**Estimated Breaking Points:**

| Scale | Nodes | Edges | Vectors | Status |
|-------|-------|-------|---------|--------|
| Small | 100K | 500K | 100K | Comfortable |
| Medium | 1M | 5M | 1M | Manageable |
| Large | 5M | 25M | 5M | **Strained** |
| XL | 10M | 50M | 10M | **Failing** |

**Single-Node Ceiling:** ~5M nodes with embeddings (Year 3-4 projection)

**Migration Path:**
1. **Immediate:** Read replicas via SQLite replication
2. **Short-term:** Shard by entity type (genes, SNPs, publications)
3. **Long-term:** Move to Kuzu or distributed system

---

### Kuzu Scaling Limits

**Architecture:** Embedded C++ with columnar storage

| Constraint | Limit | Impact |
|------------|-------|--------|
| Single-node memory | Server RAM | All hot data |
| Disk storage | Filesystem limits | Columnar compression helps |
| Query parallelism | CPU cores | Scales with hardware |
| Write concurrency | Single writer | Transaction serialization |

**Estimated Breaking Points:**

| Scale | Nodes | Edges | Status |
|-------|-------|-------|--------|
| Small | 100K | 500K | Easy |
| Medium | 1M | 5M | Comfortable |
| Large | 10M | 50M | Manageable |
| XL | 50M | 250M | **Strained** |
| XXL | 100M | 500M | **Failing** |

**Single-Node Ceiling:** ~50M nodes (Year 5+ projection, well beyond needs)

**Migration Path:**
1. **Immediate:** Vertical scaling (more RAM, faster SSD)
2. **Short-term:** External vector search (Qdrant/Weaviate)
3. **Long-term:** Consider distributed graph (Neptune, TigerGraph)

---

### Neo4j Scaling Limits

**Architecture:** JVM-based, native graph storage

| Constraint | Limit | Impact |
|------------|-------|--------|
| JVM heap | ~31 GB practical | Index structures |
| Page cache | Server RAM - heap | Property access speed |
| Single-writer | 1 thread | Write throughput |
| Causal clustering | 3+ nodes | Read scaling |

**Estimated Breaking Points (Community):**

| Scale | Nodes | Edges | Status |
|-------|-------|-------|--------|
| Small | 100K | 500K | Easy |
| Medium | 1M | 5M | Comfortable |
| Large | 10M | 50M | Manageable |
| XL | 50M | 250M | **Strained** |

**Single-Node Ceiling (Community):** ~50M nodes

**Enterprise Scaling:**
- Causal clustering for read scaling
- Fabric for federation
- Can reach 100B+ nodes with proper architecture

**Migration Path:**
1. **Immediate:** Enterprise license + read replicas
2. **Short-term:** Fabric for sharding
3. **Long-term:** Already distributed-ready

---

### SurrealDB Scaling Limits

**Architecture:** Rust-based, multi-model

| Constraint | Limit | Impact |
|------------|-------|--------|
| Single-node storage | System RAM + disk | Depends on data model |
| Concurrent connections | ~10,000 | Connection pooling helps |
| Complex queries | Memory-bound | Multi-hop expensive |
| Vector search | Newer feature | Less optimized |

**Estimated Breaking Points:**

| Scale | Nodes | Edges | Status |
|-------|-------|-------|--------|
| Small | 100K | 500K | Easy |
| Medium | 1M | 5M | Comfortable |
| Large | 5M | 25M | **Strained** |
| XL | 10M | 50M | **Failing** |

**Single-Node Ceiling:** ~5M nodes (similar to RuVector)

**Migration Path:**
1. **Immediate:** SurrealDB clustering (beta)
2. **Short-term:** Horizontal scaling via TiKV backend
3. **Long-term:** Fully distributed multi-region

---

### Scaling Ceiling Summary

| Database | Single-Node Ceiling | Year Reached* | Migration Complexity |
|----------|---------------------|---------------|---------------------|
| RuVector | ~5M nodes | Year 3-4 | Medium |
| Kuzu | ~50M nodes | Year 5+ | Low |
| Neo4j Community | ~50M nodes | Year 5+ | Low (upgrade path) |
| Neo4j Enterprise | 100B+ nodes | Never | N/A |
| SurrealDB | ~5M nodes | Year 3-4 | Low |

*Based on 50% annual data growth from 1.7M baseline

---

### Recommended Scaling Strategy

**For Gene-Disease Knowledge Graph:**

```
Year 1-2: Single-node (any option works)
    └── 1.7M → 3.8M nodes
    └── All options comfortable

Year 3-4: Optimization phase
    └── 3.8M → 8.5M nodes
    └── RuVector/SurrealDB: Consider migration
    └── Kuzu/Neo4j: Vertical scaling sufficient

Year 5+: Distributed consideration
    └── 8.5M+ nodes
    └── All options except Neo4j Enterprise need planning
```

**Recommendation:** Choose Kuzu or Neo4j Community for longest single-node runway, with clear upgrade path when needed.

---

## 5. Benchmark Methodology

### Pre-Commitment Testing Strategy

Before committing to any database, run these benchmarks with representative data.

---

### Test Dataset Specifications

**Recommended Test Scale:** 10% of production (~170K nodes)

```
Test Dataset Composition:
├── Genes: 2,000 nodes
├── SNPs: 120,000 nodes
├── Publications: 50,000 nodes
├── Pathways: 200 nodes
├── Conditions: 500 nodes
├── Herbs: 1,000 nodes
├── Nutrients: 500 nodes
└── Total: ~174,200 nodes

Edges:
├── Gene → SNP: 120,000
├── SNP → Condition: 240,000
├── Gene → Pathway: 8,000
├── Herb → Condition: 5,000
├── Publication → Entity: 150,000
└── Total: ~523,000 edges

Vectors:
├── All nodes: 174,200 × 384-dim
└── Total: ~268 MB vectors
```

---

### Core Benchmark Queries

#### Benchmark 1: Simple Lookup (Baseline)

```cypher
// Find gene by symbol
MATCH (g:Gene {symbol: 'MTHFR'})
RETURN g
```

**Target:** < 5ms
**Measure:** p50, p95, p99 latency

#### Benchmark 2: Single-Hop Traversal

```cypher
// Find all SNPs for a gene
MATCH (g:Gene {symbol: 'MTHFR'})-[:HAS_VARIANT]->(snp:SNP)
RETURN snp
```

**Target:** < 20ms
**Measure:** p50, p95, p99 latency

#### Benchmark 3: Multi-Hop Traversal (3-hop)

```cypher
// Find conditions related to gene through SNPs
MATCH (g:Gene {symbol: 'MTHFR'})
    -[:HAS_VARIANT]->(snp:SNP)
    -[:ASSOCIATED_WITH]->(c:Condition)
RETURN DISTINCT c.name, count(snp) as snp_count
ORDER BY snp_count DESC
```

**Target:** < 100ms
**Measure:** p50, p95, p99 latency, result count

#### Benchmark 4: Deep Traversal (6-hop)

```cypher
// Full pathway analysis
MATCH path = (g:Gene {symbol: 'MTHFR'})
    -[:HAS_VARIANT]->(snp:SNP)
    -[:ASSOCIATED_WITH]->(c:Condition)
    -[:TREATED_BY]->(h:Herb)
    -[:AFFECTS]->(p:Pathway)
    -[:INVOLVES]->(g2:Gene)
WHERE g2 <> g
RETURN path
LIMIT 100
```

**Target:** < 500ms
**Measure:** p50, p95, p99 latency, paths found

#### Benchmark 5: Vector Similarity Search

```javascript
// Find similar publications
const query_embedding = await embedText("MTHFR folate metabolism");
const results = await vectorSearch(query_embedding, {
  topK: 100,
  threshold: 0.3
});
```

**Target:** < 30ms for top-100
**Measure:** p50, p95, p99 latency, recall@100

#### Benchmark 6: Vector + Graph Combined

```javascript
// Semantic search with graph filter
const similar_pubs = await vectorSearch(query, { topK: 1000 });
const filtered = await graphQuery(`
  MATCH (pub:Publication)-[:MENTIONS]->(g:Gene)
  WHERE pub.id IN $ids AND g.symbol IN ['MTHFR', 'MTR']
  RETURN pub, g
`, { ids: similar_pubs.map(p => p.id) });
```

**Target:** < 100ms combined
**Measure:** Vector time, graph time, total time

#### Benchmark 7: Aggregation Query

```cypher
// Gene influence analysis
MATCH (g:Gene)-[r]-()
RETURN g.symbol, labels(g), count(r) as degree
ORDER BY degree DESC
LIMIT 100
```

**Target:** < 200ms
**Measure:** p50, p95, p99 latency

#### Benchmark 8: Batch Write

```javascript
// Insert 10,000 nodes with embeddings
await batchInsert(nodes_10k);
```

**Target:** < 60 seconds for 10K nodes
**Measure:** Throughput (nodes/second)

#### Benchmark 9: Concurrent Read Under Write

```javascript
// Simulate production load
// 90% reads, 10% writes
Promise.all([
  runReadQueries(queries, 90),  // 90 concurrent readers
  runWriteQueries(updates, 10)   // 10 concurrent writers
]);
```

**Target:** Read latency increase < 50%
**Measure:** Read p95 with/without writes

#### Benchmark 10: Cold Start Recovery

```bash
# Stop database, clear caches, restart, query immediately
systemctl stop database
sync && echo 3 > /proc/sys/vm/drop_caches
systemctl start database
time run_benchmark_query
```

**Target:** < 30 seconds to operational
**Measure:** Time to first successful query

---

### Benchmark Execution Plan

```
Week 1: Setup
├── Day 1-2: Generate test dataset
├── Day 3: Deploy RuVector test instance
├── Day 4: Deploy Kuzu test instance
├── Day 5: Deploy Neo4j test instance

Week 2: Benchmarking
├── Day 1-2: Run all benchmarks on RuVector
├── Day 3-4: Run all benchmarks on Kuzu
├── Day 5: Run all benchmarks on Neo4j

Week 3: Analysis
├── Day 1-2: Compile results
├── Day 3: Write comparison report
├── Day 4-5: Make final recommendation
```

---

### Benchmark Environment

**Recommended Test Infrastructure:**

| Component | Specification |
|-----------|---------------|
| Instance | AWS r6g.xlarge (4 vCPU, 32 GB) |
| Storage | 200 GB gp3 SSD |
| Network | Isolated VPC |
| OS | Ubuntu 22.04 LTS |
| Node.js | v20 LTS |

**Benchmark Tools:**
- k6 or artillery for load testing
- Custom scripts for database-specific queries
- Prometheus/Grafana for metrics collection

---

### Success Criteria

| Metric | Minimum | Target | Stretch |
|--------|---------|--------|---------|
| Simple lookup p95 | < 10ms | < 5ms | < 2ms |
| 3-hop traversal p95 | < 200ms | < 100ms | < 50ms |
| 6-hop traversal p95 | < 1s | < 500ms | < 200ms |
| Vector search p95 | < 50ms | < 30ms | < 10ms |
| Combined query p95 | < 200ms | < 100ms | < 50ms |
| Batch write throughput | 500 n/s | 1000 n/s | 2000 n/s |
| Cold start | < 60s | < 30s | < 10s |

---

## Summary and Recommendations

### Key Findings

1. **Memory:** All options fit within 32 GB RAM for Year 1-2
2. **Cost:** Embedded options (RuVector, Kuzu) are 3-5x cheaper than Neo4j
3. **Performance:** Neo4j leads for deep traversals, Kuzu for bulk operations
4. **Ceiling:** Kuzu and Neo4j have 10x longer runway than RuVector/SurrealDB

### Scalability Recommendation

For the Gene-Disease Knowledge Graph with target scale of ~1.7M nodes growing to ~8M over 5 years:

**Primary Recommendation: Kuzu**
- Best cost/performance ratio
- Longest single-node runway
- Clean integration path with external vector search
- Modern codebase, actively developed

**Alternative: RuVector (if TypeScript ecosystem valued)**
- Unified stack with claude-flow
- Good enough for Year 1-3
- Requires migration planning for Year 4+

**Premium Option: Neo4j Enterprise**
- If budget allows (~$50K+/year)
- Best graph algorithms (GDS library)
- Enterprise support and clustering

### Next Steps

1. Run benchmark suite on 10% dataset
2. Validate performance assumptions
3. Test vector search integration patterns
4. Document migration paths
5. Make final selection

---

**Document Status:** Complete
**Last Updated:** 2026-01-21
**Author:** Scalability & Performance Analyst (Hive-Mind)

