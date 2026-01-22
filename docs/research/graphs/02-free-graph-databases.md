# FREE Graph Database Alternatives to RuVector

> Research conducted: 2026-01-21
> Focus: Truly free/open source graph databases (Apache 2.0 / MIT / MPL 2.0)
> Excludes: Commercial licenses, BSL, SSPL, and other restrictive licenses

## Overview

This document analyzes free and open-source graph databases as potential alternatives or complements to RuVector for the Gene platform. We specifically exclude databases with commercial or restrictive licenses (Neo4j, TigerGraph, Memgraph, ArangoDB, SurrealDB, FalkorDB).

---

## Summary Table

| Database | License | Query Language | Hypergraph | Vector Search | Rust Client | Best For |
|----------|---------|----------------|------------|---------------|-------------|----------|
| Apache AGE | Apache 2.0 | Cypher + SQL | NO | NO | YES | PostgreSQL users needing graphs |
| TypeDB | MPL 2.0 | TypeQL | YES | NO | YES | Knowledge graphs, semantic modeling |
| CozoDB | MPL 2.0 | Datalog | NO | YES (HNSW) | Native Rust | AI/ML, vector+graph hybrid |
| NebulaGraph | Apache 2.0 | nGQL/Cypher | NO | NO | YES | Large-scale distributed graphs |
| Oxigraph | Apache 2.0/MIT | SPARQL | NO | NO | Native Rust | RDF/semantic web |
| IndraDB | MPL 2.0 | Custom API | NO | NO | Native Rust | Embedded graph apps |
| JanusGraph | Apache 2.0 | Gremlin | NO | NO | YES (gremlin-rs) | Enterprise distributed graphs |
| agdb | Apache 2.0 | Native API | NO | NO | Native Rust | Embedded ACID graph |

---

## 1. Apache AGE

**License:** Apache 2.0 (Truly Free)
**Query Language:** openCypher + SQL hybrid
**Website:** https://age.apache.org/

### Overview

Apache AGE (A Graph Extension) is a PostgreSQL extension that provides graph database functionality. It enables users to leverage graph capabilities on top of existing PostgreSQL infrastructure.

### Key Features

- **PostgreSQL Integration:** Supports PostgreSQL 11-17
- **Hybrid Queries:** Cypher queries can be embedded within SQL using a `cypher()` function
- **Property Graph Model:** Supports vertices, edges, and properties
- **Custom Data Type:** Uses `agtype` (superset of JSON) for all graph data
- **Full SQL Compatibility:** Use standard ANSI SQL alongside openCypher

### Cypher Support

Supported operations:
- MATCH, WITH, RETURN, ORDER BY, SKIP, LIMIT
- CREATE, DELETE, SET, REMOVE, MERGE
- Predicate, Scalar, List, Numeric, String, Aggregation Functions

### Limitations

- **No Hypergraph Support:** Binary edges only (not n-ary)
- **No Vector Search:** Must use pgvector separately
- **Cypher in Expressions:** `ERROR: cypher(...) in expressions is not supported` - must use subqueries
- **Null Handling:** null is not equal to null (Cypher semantics)
- **Performance:** Graph traversals may not be as optimized as native graph DBs

### Rust Client

- **Crate:** `apache_age` (v0.6.4)
- **GitHub:** https://github.com/Dzordzu/rust-apache-age
- **Based on:** postgres and tokio-postgres crates
- **Features:** Sync and async clients

### Best Use Cases

- Organizations already using PostgreSQL
- Hybrid relational + graph workloads
- Migration path from relational to graph

### Sources

- [Apache AGE Official](https://age.apache.org/)
- [GitHub - apache/age](https://github.com/apache/age)
- [Rust Driver](https://crates.io/crates/apache_age)

---

## 2. TypeDB

**License:** MPL 2.0 (Truly Free - Community Edition)
**Query Language:** TypeQL
**Website:** https://typedb.com/

### Overview

TypeDB is a next-generation database with a hypergraph data model and polymorphic type system. Version 3.0 (December 2024) was rewritten in Rust.

### Key Features

- **Hypergraph Model:** Relations connect any number of entities (n-ary, not binary)
- **Polymorphic Types:** Entities, relations, and attributes are all subtypable
- **TypeQL:** Declarative, functional, strongly-typed query language
- **Recursive Functions:** Datalog-like functions with recursion support
- **Inference Engine:** Built-in reasoning capabilities

### Hypergraph Support - YES

TypeDB's data structure is a true hypergraph:
- Relations can connect any number of entities or attributes
- Relations can be nested (a hyperedge is also a vertex)
- Entities and relations can own multiple attributes of same type
- Maintains true context without artificial joins

### Schema and Type System

- Polymorphic Entity-Relation-Attribute (PERA) model
- Inheritance and interfaces
- Subtypes inherit behaviors from supertypes
- No normalization or reification required

### Limitations

- **No Vector Search:** Must use external vector store
- **Learning Curve:** TypeQL is different from SQL/Cypher
- **Performance:** May be slower than simpler graph DBs for basic traversals
- **Cloud Features:** Some features only in commercial TypeDB Cloud

### Rust Client

- **Crate:** `typedb-driver`
- **GitHub:** https://github.com/typedb/typedb-driver-rust
- **Features:** Fully async API, sync via feature flag, token-based auth
- **Min Rust:** 1.68.0

### Best Use Cases

- Knowledge graphs with complex relationships
- Semantic data modeling
- Applications requiring inference/reasoning
- Biomedical and scientific data
- RAG systems for AI

### Sources

- [TypeDB Official](https://typedb.com/)
- [TypeDB GitHub](https://github.com/typedb/typedb)
- [Hypergraph Modeling](https://medium.com/vaticle/modelling-data-with-hypergraphs-edff1e12edf0)
- [Rust Driver](https://crates.io/crates/typedb-driver)

---

## 3. CozoDB

**License:** MPL 2.0 (Truly Free)
**Query Language:** Datalog
**Website:** https://www.cozodb.org/

### Overview

CozoDB is a transactional, relational-graph-vector database that uses Datalog for queries. Marketed as "The hippocampus for AI!"

### Key Features

- **Hybrid Database:** Relational + Graph + Vector in one
- **Datalog Queries:** Powerful recursive queries, highly composable
- **HNSW Vector Search:** Integrated with Datalog queries
- **Time Travel:** Query historical data states
- **Written in Rust:** High performance, no GC pauses

### Vector Search Support - YES (HNSW)

HNSW indices integrate seamlessly with Datalog:
- Create HNSW indices on relations containing vectors
- Multiple indices per relation with filters
- Use vector search in recursive Datalog
- Distance metrics: L2 (default), Cosine, IP
- MVCC protection for concurrent writes
- Memory-efficient: data resides on disk

### Additional Features (v0.7+)

- MinHash-LSH for near-duplicate search
- Full-text search
- JSON value support

### Limitations

- **No Hypergraph Support:** Standard graph model
- **Datalog Learning Curve:** Different from SQL/Cypher
- **Smaller Community:** Less ecosystem than major graph DBs
- **Documentation:** Less comprehensive than established DBs

### Rust Client

- **Native Rust:** CozoDB is written in Rust, use directly as library
- **Crate:** `cozo` on crates.io
- **Embeddable:** Can embed directly in Rust applications

### Best Use Cases

- AI/ML applications needing vector + graph
- Complex recursive graph queries
- LLM memory systems ("hippocampus")
- Analytics requiring time travel
- Embedding directly in Rust applications

### Sources

- [CozoDB Official](https://www.cozodb.org/)
- [GitHub - cozodb/cozo](https://github.com/cozodb/cozo)
- [Vector Search Docs](https://docs.cozodb.org/en/latest/vector.html)

---

## 4. NebulaGraph

**License:** Apache 2.0 (Truly Free)
**Query Language:** nGQL (SQL-like), OpenCypher compatible
**Website:** https://www.nebula-graph.io/

### Overview

NebulaGraph is a distributed graph database designed for super large-scale graphs with milliseconds of latency. Claims to be the only open-source graph DB that can handle trillions of edges.

### Key Features

- **Distributed Architecture:** Shared-nothing, linear scalability
- **Three Services:** Graph, Storage (Raft), Meta
- **High Performance:** Millisecond latency, high throughput
- **nGQL:** Easy SQL-like query language
- **OpenCypher/GQL:** Compatibility for migration

### Architecture

- **Graph Service:** Stateless query processing (nebula-graphd)
- **Storage Service:** Distributed Raft-based storage (nebula-storaged)
- **Meta Service:** Schema, user accounts, job management

### Limitations

- **No Hypergraph Support:** Property graph model only
- **No Vector Search:** Must use external vector store
- **Complexity:** Requires multiple services for production
- **Resource Heavy:** Needs significant infrastructure for distributed mode

### Rust Client

- **Official:** `nebula-rust` by vesoft-inc
- **Community:** `rust-nebula` by nebula-contrib (V3 protocol)
- **Crates:** `nebula-client`, `nebula-graph-client`
- **Features:** Graph, Meta, and Storage clients

### Best Use Cases

- Social networks at scale
- Recommendation systems
- Knowledge graphs (large scale)
- Fraud detection
- Real-time graph analytics

### Production Users

Snapchat, Binance, Akulaku, GrowingIO, IRL

### Sources

- [NebulaGraph Official](https://www.nebula-graph.io/)
- [GitHub - vesoft-inc/nebula](https://github.com/vesoft-inc/nebula)
- [Rust Client](https://github.com/nebula-contrib/rust-nebula)

---

## 5. Oxigraph

**License:** Apache 2.0 / MIT (Dual License - Truly Free)
**Query Language:** SPARQL 1.1 (+ RDF 1.2 support)
**Website:** https://github.com/oxigraph/oxigraph

### Overview

Oxigraph is a graph database implementing the SPARQL standard, written in Rust. Designed for RDF (Resource Description Framework) data and semantic web applications.

### Key Features

- **SPARQL 1.1 Full Support:** Query, Update, Federated Query
- **RDF 1.2 Support:** Behind `rdf-12` feature flag
- **RDF Formats:** Turtle, TriG, N-Triples, N-Quads, RDF/XML
- **Storage Options:** RocksDB (on-disk) or in-memory
- **WASM Support:** Can compile to WebAssembly
- **CLI Tool:** Standalone server with YASGUI web UI

### Recent Updates (2025)

- RDF 1.2 support (replacing RDF-star)
- Custom aggregate functions in Rust/Python API
- CancellationToken for query timeout
- Parallel parsing for NTriples/NQuads

### Limitations

- **No Hypergraph Support:** RDF triple model (subject-predicate-object)
- **No Vector Search:** Pure RDF/SPARQL focus
- **RDF-Specific:** Not a general property graph database
- **Query Optimization:** Not fully optimized yet (per documentation)

### Rust Client

- **Native Rust:** Oxigraph is written in Rust
- **Crate:** `oxigraph` on crates.io
- **Modular:** `oxrdf`, `oxrdfio`, `oxttl` as standalone crates
- **Embeddable:** Use as library in Rust applications

### Best Use Cases

- Semantic web applications
- RDF data management
- SPARQL endpoints
- Linked data applications
- Ontology-based systems

### Sources

- [GitHub - oxigraph/oxigraph](https://github.com/oxigraph/oxigraph)
- [Rust Docs](https://docs.rs/oxigraph)
- [Crates.io](https://crates.io/crates/oxigraph)

---

## 6. IndraDB

**License:** MPL 2.0 (Truly Free)
**Query Language:** Custom API (no query language)
**Website:** https://indradb.github.io/

### Overview

IndraDB is a graph database written in Rust, originally inspired by Facebook's TAO datastore. Emphasizes simplicity and is designed for graphs too large for full processing.

### Key Features

- **Directed Typed Graphs:** Support for multi-hop queries
- **JSON Properties:** On vertices and edges
- **Pluggable Datastores:** In-memory, PostgreSQL, sled
- **gRPC Server:** Cross-language support
- **Embeddable:** Use as Rust library directly

### TAO Heritage

Originally TAO-like design focusing on:
- Simplicity of implementation
- Simple query semantics
- Handling graphs too large for full processing

Note: Now has richer features, no longer strictly TAO-like.

### Limitations

- **No Hypergraph Support:** Standard directed graph model
- **No Vector Search:** No built-in vector capabilities
- **No Query Language:** API-based only (no Cypher/SQL/SPARQL)
- **Limited Features:** Simpler than feature-rich graph DBs
- **Smaller Community:** Less documentation and ecosystem

### Rust Client

- **Native Rust:** IndraDB is written in Rust
- **Crate:** `indradb` (v4.0.0)
- **Server Crate:** `indradb-proto` for gRPC
- **Bindings:** Python and Rust official bindings

### Best Use Cases

- Embedding graph DB in Rust applications
- Simple graph storage needs
- Applications needing TAO-style semantics
- Microservices with gRPC

### Sources

- [IndraDB Official](https://indradb.github.io/)
- [GitHub - indradb/indradb](https://github.com/indradb/indradb)
- [Crates.io](https://crates.io/crates/indradb)

---

## 7. JanusGraph

**License:** Apache 2.0 (Truly Free)
**Query Language:** Gremlin (Apache TinkerPop)
**Website:** https://janusgraph.org/

### Overview

JanusGraph is an open-source, distributed graph database under The Linux Foundation. Designed for graphs requiring storage beyond single-machine capacity.

### Key Features

- **Gremlin Query Language:** Apache TinkerPop standard
- **Multiple Backends:** Cassandra, HBase, Bigtable, BerkeleyDB, ScyllaDB
- **Full-Text Search:** Elasticsearch, Solr, Lucene integration
- **Big Data Integration:** Spark, Giraph, Hadoop
- **No Vendor Lock-in:** Standard Gremlin enables migration

### Architecture

- Uses pluggable storage backends
- Separates storage from compute
- Supports distributed transactions
- Multi-datacenter replication (with Cassandra)

### Limitations

- **No Hypergraph Support:** Property graph model only
- **No Vector Search:** Must use external systems
- **Complexity:** Requires storage backend setup
- **Java-Based:** Written in Java (not Rust-native)
- **Operational Overhead:** Multiple components to manage

### Rust Client

- **gremlin-rs:** Experimental Rust client for TinkerPop
- **Crate:** `gremlin-client`
- **Features:** Sync and async, GLV traversal support
- **Derive Macros:** `FromGValue`, `FromGMap` for deserialization

### Best Use Cases

- Enterprise-scale graph analytics
- Systems already using Cassandra/HBase
- Big data graph processing
- Applications requiring Gremlin standard

### Production Support

Supported by IBM, Google, Hortonworks

### Sources

- [JanusGraph Official](https://janusgraph.org/)
- [GitHub - JanusGraph/janusgraph](https://github.com/JanusGraph/janusgraph)
- [gremlin-rs](https://github.com/wolf4ood/gremlin-rs)

---

## 8. agdb (Agnesoft Graph Database)

**License:** Apache 2.0 (Truly Free)
**Query Language:** Native API (no query language)
**Website:** https://agnesoft.com/en-US/whyagdb

### Overview

agdb is an application-native, embedded graph database written in safe Rust. Focuses on minimal hardware footprint with full ACID compliance.

### Key Features

- **Full ACID:** Memory-mapped file with write-ahead log
- **No Query Language:** Queries in same language as application
- **Schema-less:** Typed but flexible data store
- **Constant Time Operations:** Lookups/traversals regardless of size
- **Cluster Mode:** Custom consensus protocol for distributed use

### Storage Options

- **Default (Db):** Memory-mapped, full ACID, faster
- **File-Based (DbFile):** Minimum memory footprint, ACID with WAL
- Files are interchangeable between modes

### Limitations

- **No Hypergraph Support:** Standard graph model
- **No Vector Search:** No built-in vector capabilities
- **No Query Language:** Must use Rust API directly
- **Smaller Community:** Newer, less established
- **Rust-Only:** No cross-language support (except via bindings)

### Rust Client

- **Native Rust:** agdb is written in Rust
- **Crate:** `agdb` (v0.11.2)
- **Features:** Safe Rust, no unsafe code

### Best Use Cases

- Embedded graph applications
- Rust-native development
- Applications needing ACID guarantees
- In-memory caching with persistence
- Microservices with embedded state

### Sources

- [agdb Docs](https://docs.rs/agdb/latest/agdb/)
- [GitHub - agnesoft/agdb](https://github.com/agnesoft/agdb)
- [Why agdb?](https://agnesoft.com/en-US/whyagdb)

---

## Recommendations for Gene Platform

### For Vector + Graph (RuVector Alternative)

**CozoDB** is the strongest candidate:
- Native HNSW vector search integrated with Datalog
- True hybrid: relational + graph + vector
- Written in Rust, embeddable
- MPL 2.0 license (truly free)

### For Hypergraph/Knowledge Graph

**TypeDB** is the only option with true hypergraph support:
- N-ary relations (not just binary edges)
- Polymorphic type system with inference
- TypeQL for complex semantic queries
- MPL 2.0 license (Community Edition)

### For PostgreSQL Integration

**Apache AGE** for existing PostgreSQL infrastructure:
- Leverage existing PostgreSQL investment
- Hybrid SQL + Cypher queries
- Apache 2.0 license

### For Large-Scale Distributed

**NebulaGraph** or **JanusGraph**:
- NebulaGraph: Faster, nGQL, trillions of edges
- JanusGraph: Gremlin standard, flexible backends
- Both Apache 2.0

### For Embedded Rust Applications

**CozoDB**, **Oxigraph**, **IndraDB**, or **agdb**:
- All native Rust
- CozoDB: Best for vector + graph hybrid
- Oxigraph: Best for RDF/SPARQL
- IndraDB: Simple property graph
- agdb: Minimal footprint, ACID

---

## Excluded Databases (Commercial/Restrictive Licenses)

The following were explicitly excluded from this analysis:

| Database | License | Reason for Exclusion |
|----------|---------|---------------------|
| Neo4j | Commercial (AGPL for Community) | AGPL restrictions on commercial use |
| TigerGraph | Commercial | Proprietary license |
| Memgraph | BSL 1.1 | Business Source License restrictions |
| ArangoDB | BSL 1.1 | Business Source License restrictions |
| SurrealDB | BSL 1.1 | Business Source License restrictions |
| FalkorDB | SSPLv1 | Server Side Public License restrictions |

---

## Conclusion

For RuVector alternatives in the Gene platform:

1. **Best Overall:** CozoDB - only free option with both graph AND vector search
2. **Best for Knowledge Graphs:** TypeDB - only free option with true hypergraph support
3. **Best for Scale:** NebulaGraph - proven at trillion-edge scale
4. **Best for PostgreSQL Users:** Apache AGE - extend existing infrastructure

Consider using **CozoDB** as primary if vector search is critical, potentially combined with **TypeDB** for complex semantic modeling requirements.
