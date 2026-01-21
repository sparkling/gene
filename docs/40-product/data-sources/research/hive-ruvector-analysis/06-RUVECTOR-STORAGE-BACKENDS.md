# RuVector Storage Backends Analysis

**Research Date:** 2025-01-21
**Analyst:** Hive-Mind RuVector Storage Backend Specialist
**Status:** COMPLETE

---

## Executive Summary

RuVector is a distributed vector database built in Rust that offers multiple storage backend options to suit different deployment scenarios. This document provides comprehensive analysis of each storage backend available in RuVector.

**Key Finding:** RuVector uses a **tiered storage architecture** that automatically manages data placement based on access patterns - hot data stays in fast storage, cold data moves to denser storage.

---

## 1. Overview of RuVector Storage Architecture

RuVector's storage is designed with a hierarchical approach similar to a computer's memory hierarchy:

- **Hot Tier:** Frequently accessed vectors stay in fast cache/memory
- **Warm Tier:** Moderate access data in efficient storage
- **Cold Tier:** Rarely accessed data compressed in dense storage

This tiering happens **automatically with no configuration needed** - RuVector tracks access patterns and promotes/demotes vectors between tiers.

### Core Storage Dependencies (from Cargo.toml)

Based on research, RuVector's storage layer uses:
- **redb** - Embedded key-value store for persistent storage
- **memmap2** - Memory-mapped files for efficient I/O
- **PostgreSQL** (via ruvector-postgres extension) - For SQL-based deployments

---

## 2. Storage Backend: redb (Embedded Key-Value Store)

### What is redb?

[redb](https://github.com/cberner/redb) is a simple, portable, high-performance, ACID, embedded key-value store written in pure Rust. It's designed to be a modern alternative to LMDB with better safety guarantees.

### How RuVector Uses redb

RuVector uses redb as its **primary persistent storage backend** for:
- Vector data storage
- Index metadata persistence
- Graph relationship storage
- Configuration and state

### Configuration Options

```rust
use ruvector_core::{VectorDB, DbOptions};

let mut options = DbOptions::default();
options.dimensions = 128;
options.storage_path = "./vectors.db".to_string();  // redb database file
let db = VectorDB::new(options)?;
```

### Performance Characteristics

| Metric | Value |
|--------|-------|
| Write Performance | High throughput, ACID compliant |
| Read Latency | Sub-millisecond for cached data |
| Durability | Full ACID transactions |
| Memory Usage | Efficient with memory-mapping |
| Max Database Size | Limited by disk space |

### When to Use redb

- **Embedded/Edge deployments** - Single-node applications
- **CLI tools** - Command-line vector search
- **Development/Testing** - Local development environments
- **Serverless** - Stateless functions with persistent state
- **IoT/ESP32** - Embedded systems (RuVector supports ESP32)

### Advantages
- Zero external dependencies
- Pure Rust, memory-safe
- ACID transactions
- Portable across platforms
- No server process needed

### Limitations
- Single-process access (no concurrent writers from multiple processes)
- Not suitable for high-concurrency multi-process scenarios

---

## 3. Storage Backend: memmap2 (Memory-Mapped Files)

### What is memmap2?

[memmap2](https://crates.io/crates/memmap2) is a Rust library for memory-mapped file I/O, allowing files to be accessed as if they were in memory.

### How RuVector Uses memmap2

RuVector uses memory-mapped files for:
- **Large vector arrays** - Efficient access without loading entire datasets
- **HNSW graph structures** - Fast graph traversal
- **Index data** - Quick lookup without serialization overhead

### Architecture

```
┌─────────────────────────────────────────┐
│           Virtual Memory Space           │
├─────────────────────────────────────────┤
│  Memory-Mapped Vector Data (memmap2)    │
│  ┌─────┬─────┬─────┬─────┬─────┐       │
│  │ V1  │ V2  │ V3  │ ... │ Vn  │       │
│  └─────┴─────┴─────┴─────┴─────┘       │
├─────────────────────────────────────────┤
│  Memory-Mapped HNSW Graph               │
│  (Node connections, layer data)         │
├─────────────────────────────────────────┤
│  redb Metadata Store                    │
└─────────────────────────────────────────┘
           │
           ▼
    ┌─────────────┐
    │  Disk Files │
    └─────────────┘
```

### Performance Characteristics

| Metric | Benefit |
|--------|---------|
| Startup Time | Near-instant (lazy loading) |
| Memory Efficiency | Only active pages loaded |
| Large Dataset Support | Handles datasets larger than RAM |
| OS Page Cache | Leverages kernel's caching |

### Configuration

Memory mapping is typically automatic in RuVector but can be influenced by:
- File location (SSD vs HDD significantly impacts performance)
- Available system memory for page cache
- Access patterns (sequential vs random)

### When to Use memmap2 Backend

- **Large vector collections** - Millions to billions of vectors
- **Memory-constrained environments** - Less RAM than dataset size
- **Read-heavy workloads** - Optimized for search queries
- **Fast startup required** - No loading delay

---

## 4. Storage Backend: PostgreSQL (ruvector-postgres)

### Overview

**ruvector-postgres** is a pgvector-compatible PostgreSQL extension providing:
- Drop-in replacement for pgvector
- 77+ SQL functions for vector operations
- Full SIMD acceleration (AVX-512/AVX2/NEON)
- Self-learning capabilities

### Version: 0.2.0

### Installation Options

```bash
# Docker (recommended)
docker run -d -e POSTGRES_PASSWORD=secret -p 5432:5432 ruvector/postgres:latest

# From source
cargo install cargo-pgrx --version "0.12.9" --locked
cargo pgrx install --release

# CLI tool for management
npm install -g @ruvector/postgres-cli
ruvector-pg install
ruvector-pg vector create table --dim 1536 --index hnsw
```

### Key Features

| Feature | Description |
|---------|-------------|
| **77+ SQL Functions** | Comprehensive vector operations |
| **SIMD Acceleration** | ~2x faster than standard AVX2 |
| **Index Types** | HNSW and IVFFlat |
| **39 Attention Mechanisms** | Neural attention patterns |
| **GNN Layers** | Graph Neural Network support |
| **Hyperbolic Embeddings** | Poincare + Lorentz models |
| **Sparse Vectors/BM25** | Hybrid search support |
| **SPARQL 1.1** | 50+ RDF functions |
| **Local Embeddings** | 6 fastembed models built-in |

### SQL Function Categories (77+ Total)

Based on the feature set, the SQL functions are organized into these categories:

| Category | Approx. Count | Examples |
|----------|---------------|----------|
| **Vector Operations** | 15+ | Distance calculations, similarity search |
| **Index Management** | 10+ | HNSW/IVFFlat create, configure, optimize |
| **Attention Mechanisms** | 39 | Self-attention, cross-attention, multi-head |
| **GNN Layers** | 5+ | Graph convolutions, message passing |
| **Hyperbolic Embeddings** | 5+ | Poincare distance, Lorentz transforms |
| **SPARQL/RDF** | 50+ | Triple patterns, graph queries |
| **Embedding Generation** | 6+ | fastembed model inference |

### PostgreSQL vs Embedded Comparison

| Aspect | redb (Embedded) | PostgreSQL (ruvector-postgres) |
|--------|-----------------|-------------------------------|
| **Setup Complexity** | Minimal | Requires PostgreSQL server |
| **Scalability** | Single-node | Multi-node with replication |
| **SQL Support** | None | Full SQL + extensions |
| **Concurrent Access** | Limited | Excellent |
| **Transactions** | ACID | Full PostgreSQL ACID |
| **Use Case** | Edge/CLI/Embedded | Production/Enterprise |

### When to Use PostgreSQL Backend

- **Enterprise deployments** - Need SQL interface and DBA familiarity
- **Existing PostgreSQL infrastructure** - Add vector capabilities
- **High concurrency** - Multiple writers and readers
- **Complex queries** - JOIN with relational data
- **Backup/Recovery** - Use standard PostgreSQL tooling

---

## 5. Storage Backend: In-Memory (Volatile)

### Overview

RuVector supports pure in-memory storage for maximum performance when persistence is not required.

### Configuration

**Rust:**
```rust
use ruvector_core::{VectorDB, DbOptions};

let mut options = DbOptions::default();
options.dimensions = 384;
// Omit storage_path for in-memory operation
let db = VectorDB::new(options)?;
```

**JavaScript/TypeScript (npm):**
```javascript
import { VectorDB } from 'ruvector';

// In-memory only - no storagePath specified
const db = new VectorDB({
  dimensions: 384,
  metric: 'cosine'
});

// For persistence, add storagePath
const persistentDb = new VectorDB({
  dimensions: 384,
  metric: 'cosine',
  storagePath: './persistent.db'  // Enables auto-save to disk
});
```

### Performance Characteristics

| Metric | Value |
|--------|-------|
| Query Latency | <0.5ms p50 |
| Throughput | 50K+ ops/sec (native) |
| Startup Time | Instant |
| Data Durability | None (volatile) |

### When to Use In-Memory

- **Testing and development** - Fast iteration
- **Ephemeral workloads** - Lambda/serverless with external persistence
- **Caching layer** - Hot data cache in front of persistent store
- **Benchmarking** - Measure raw performance

---

## 6. Storage Backend: WebAssembly/Browser (IndexedDB)

### Overview

RuVector WASM enables browser-based vector search with **IndexedDB persistence** for offline-first applications.

### Architecture

```
┌─────────────────────────────────────────┐
│           Browser Environment            │
├─────────────────────────────────────────┤
│  RuVector WASM Module                   │
│  ┌─────────────────────────────────┐    │
│  │  VectorDB Instance (Memory)     │    │
│  │  - HNSW Index                   │    │
│  │  - Vector Data                  │    │
│  └─────────────────────────────────┘    │
│              │                           │
│              ▼                           │
│  ┌─────────────────────────────────┐    │
│  │  IndexedDB Persistence Layer    │    │
│  │  - Async save/load              │    │
│  │  - Offline-first support        │    │
│  └─────────────────────────────────┘    │
├─────────────────────────────────────────┤
│  Web Workers (Optional)                  │
│  - Background processing                 │
│  - Non-blocking queries                  │
└─────────────────────────────────────────┘
```

### Usage Example

```javascript
import init, { VectorDB } from '@ruvector/wasm';

// Initialize WASM module
await init();

// Create database
const db = new VectorDB({
  dimensions: 384,
  metric: 'cosine',
  hnswM: 16,
  hnswEfConstruction: 100,
  hnswEfSearch: 64
});

// Save to IndexedDB for persistence
await db.save_to_indexed_db();

// Load previously saved database
const loadedDb = await VectorDB.load_from_indexed_db();
```

### Performance (WASM)

| Environment | Latency | Throughput |
|-------------|---------|------------|
| Native (Rust/Node) | <0.5ms | 50K+ ops/sec |
| WASM (Browser) | 10-50ms | ~1K ops/sec |

### When to Use WASM/Browser Backend

- **Privacy-first applications** - Data stays in browser
- **Offline-first PWAs** - Works without network
- **Edge computing** - Client-side AI/ML
- **Demos and prototypes** - No backend required

---

## 7. Storage Backend: Distributed (Raft Consensus)

### Overview

RuVector supports **distributed storage** with Raft consensus for horizontal scaling and fault tolerance.

### Architecture

```
┌─────────────────────────────────────────────────────────┐
│                    Raft Cluster (5 nodes)                │
├─────────────────────────────────────────────────────────┤
│                                                          │
│  ┌─────────┐    ┌─────────┐    ┌─────────┐             │
│  │ Node 1  │◄──►│ Node 2  │◄──►│ Node 3  │             │
│  │(Leader) │    │(Follower│    │(Follower│             │
│  └────┬────┘    └────┬────┘    └────┬────┘             │
│       │              │              │                   │
│       ▼              ▼              ▼                   │
│  ┌─────────┐    ┌─────────┐    ┌─────────┐             │
│  │  redb   │    │  redb   │    │  redb   │             │
│  │ Storage │    │ Storage │    │ Storage │             │
│  └─────────┘    └─────────┘    └─────────┘             │
│                                                          │
│  ┌─────────┐    ┌─────────┐                             │
│  │ Node 4  │◄──►│ Node 5  │    Consistent Hash Ring    │
│  │(Follower│    │(Follower│    for data sharding       │
│  └─────────┘    └─────────┘                             │
│                                                          │
└─────────────────────────────────────────────────────────┘
```

### Configuration Example

```rust
use ruvector_raft::RaftNode;
use ruvector_cluster::{ClusterManager, ConsistentHashRing};
use ruvector_replication::SyncManager;

// Configure 5-node Raft cluster
let config = RaftConfig {
    node_id: 1,
    cluster_members: vec![
        "node1:9000", "node2:9000", "node3:9000",
        "node4:9000", "node5:9000"
    ],
    election_timeout_min: 150,
    election_timeout_max: 300,
    heartbeat_interval: 50,
};

let raft_node = RaftNode::new(config)?;
let cluster = ClusterManager::new(raft_node);
```

### Key Components

| Component | Crate | Purpose |
|-----------|-------|---------|
| **Raft Consensus** | `ruvector_raft` | Leader election, log replication |
| **Cluster Management** | `ruvector_cluster` | Node discovery, health checks |
| **Consistent Hashing** | `ruvector_cluster` | Data distribution, sharding |
| **Replication** | `ruvector_replication` | Sync across nodes |

### When to Use Distributed Backend

- **High availability** - Survive node failures
- **Large scale** - Billions of vectors
- **Geographic distribution** - Multi-region deployments
- **Mission critical** - Financial, healthcare applications

---

## 8. Automatic Tiered Storage

### Overview

RuVector implements **automatic tiered storage** that requires no configuration. The system intelligently manages data placement based on access patterns.

### How It Works

```
Access Frequency
      ▲
      │
  Hot │  ┌─────────────────────────────────────────┐
      │  │ Full precision, maximum compute         │
      │  │ In-memory cache, immediate access       │
      │  └─────────────────────────────────────────┘
      │
 Warm │  ┌─────────────────────────────────────────┐
      │  │ Efficient storage, balanced access      │
      │  │ Memory-mapped files, fast retrieval     │
      │  └─────────────────────────────────────────┘
      │
 Cold │  ┌─────────────────────────────────────────┐
      │  │ Compressed automatically, minimal CPU   │
      │  │ Background optimization, archival       │
      │  └─────────────────────────────────────────┘
      │
      └──────────────────────────────────────────────► Time
                    (Age of data)
```

### Automatic Behaviors

| Behavior | Description |
|----------|-------------|
| **Access Tracking** | Monitors query patterns per vector |
| **Auto-Promotion** | Frequently accessed data moves to hot tier |
| **Auto-Demotion** | Stale data compresses and moves to cold tier |
| **Background Optimization** | Historical data optimizes without intervention |
| **Precision Adaptation** | Hot paths get full precision; cold paths compress |

### Benefits

- **Zero configuration** - Works out of the box
- **Cost optimization** - Reduces storage costs automatically
- **Performance optimization** - Hot data stays fast
- **Resource efficiency** - Cold data minimizes CPU/memory

---

## 9. Storage Configuration Summary

### Cargo.toml Dependencies

```toml
[dependencies]
# Core storage
redb = { version = "^2.1", optional = true }
memmap2 = { version = "^0.9", optional = true }

# Serialization
bincode = "^2.0.0-rc.3"
rkyv = "^0.8"  # Zero-copy serialization

[features]
default = ["simd", "uuid-support"]
simd = []
uuid-support = []
persistence = ["redb", "memmap2"]
```

### Feature Flags

| Feature | Purpose |
|---------|---------|
| `simd` | SIMD-optimized distance calculations |
| `uuid-support` | UUID generation for vector IDs |
| `persistence` | Enable redb + memmap2 storage |

---

## 10. Backend Selection Guide

### Decision Matrix

| Requirement | Recommended Backend |
|-------------|---------------------|
| Single-node, embedded | redb + memmap2 |
| Browser/client-side | WASM + IndexedDB |
| Enterprise/SQL needed | PostgreSQL (ruvector-postgres) |
| High availability | Distributed (Raft) |
| Maximum performance | In-memory |
| Large datasets (>RAM) | memmap2 |
| Offline-first app | WASM + IndexedDB |
| Testing/development | In-memory |

### Comparison Table

| Backend | Persistence | Scalability | Latency | Setup Complexity |
|---------|-------------|-------------|---------|------------------|
| **In-Memory** | None | Single-node | <0.5ms | Minimal |
| **redb** | Full ACID | Single-node | <1ms | Low |
| **memmap2** | File-based | Single-node | <1ms | Low |
| **PostgreSQL** | Full | Multi-node | 1-5ms | Medium |
| **WASM/IndexedDB** | Browser local | Single-client | 10-50ms | Low |
| **Distributed (Raft)** | Replicated | Horizontal | 5-20ms | High |

---

## 11. Migration Between Backends

### Export/Import Pattern

```rust
// Export from one backend
let data = db.export_all()?;

// Import to another backend
let new_db = VectorDB::new(new_options)?;
new_db.import_all(data)?;
```

### Considerations

- **Data format** - RuVector uses bincode/rkyv for efficient serialization
- **Index rebuilding** - HNSW indexes may need rebuild on import
- **Metadata preservation** - Ensure metadata tags/filters transfer

---

## 12. Code Examples

### Initialize with redb Persistence

```rust
use ruvector_core::{VectorDB, DbOptions, DistanceMetric};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut options = DbOptions::default();
    options.dimensions = 384;
    options.storage_path = "./my_vectors.db".to_string();
    options.distance_metric = DistanceMetric::Cosine;

    // HNSW tuning
    options.hnsw_m = 32;  // Connections per node (16-64)
    options.hnsw_ef_construction = 200;  // Build quality
    options.hnsw_ef_search = 100;  // Query quality

    let db = VectorDB::new(options)?;
    Ok(())
}
```

### Memory-Mapped Large Dataset

```rust
let mut options = DbOptions::default();
options.dimensions = 768;
options.storage_path = "./large_dataset.db".to_string();
options.mmap_vectors = true;  // Enable memory mapping for large datasets

let db = VectorDB::new(options)?;
```

### Node.js with Persistence

```javascript
import { VectorDB } from 'ruvector';

const db = new VectorDB({
  dimensions: 384,
  metric: 'cosine',
  storagePath: './vectors.db',  // Enables persistence
  hnswM: 32,
  hnswEfConstruction: 200,
  hnswEfSearch: 100
});

// Data auto-persists
await db.add('vec-1', embedding, { category: 'product' });

// Check stats
const stats = db.stats();
console.log(`Vectors: ${stats.totalVectors}, Memory: ${stats.memoryUsage}`);
```

---

## Sources

- [RuVector GitHub Repository](https://github.com/ruvnet/ruvector)
- [ruvector-core on crates.io](https://crates.io/crates/ruvector-core)
- [ruvector-core documentation](https://docs.rs/ruvector-core/latest/ruvector_core/)
- [ruvector npm package](https://www.npmjs.com/package/ruvector)
- [ruvector-postgres on crates.io](https://crates.io/crates/ruvector-postgres)
- [redb - Embedded Database](https://github.com/cberner/redb)
- [memmap2 crate](https://crates.io/crates/memmap2)

---

**Research Status:** COMPLETE
**Last Updated:** 2025-01-21
