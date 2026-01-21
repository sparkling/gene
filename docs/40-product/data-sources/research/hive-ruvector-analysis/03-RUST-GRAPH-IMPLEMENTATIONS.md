# Rust Graph Implementations Analysis

> Research conducted: 2026-01-21
> Purpose: Evaluate Rust graph libraries for integration with RuVector and Gene platform

## Table of Contents
1. [petgraph - Standard Graph Library](#1-petgraph---standard-graph-library)
2. [graph crate - Alternative Implementation](#2-graph-crate---alternative-implementation)
3. [differential-dataflow - Timely Dataflow](#3-differential-dataflow---timely-dataflow)
4. [pathfinding crate - Algorithm Focus](#4-pathfinding-crate---algorithm-focus)
5. [Rust Graph Databases](#5-rust-graph-databases)
6. [RuVector's Internal Graph Usage](#6-ruvectors-internal-graph-usage)

---

## 1. petgraph - Standard Graph Library

**The de facto standard graph library for Rust**

| Attribute | Value |
|-----------|-------|
| Crates.io | [petgraph](https://crates.io/crates/petgraph) |
| Documentation | [docs.rs/petgraph](https://docs.rs/petgraph/) |
| GitHub | [petgraph/petgraph](https://github.com/petgraph/petgraph) |
| Current Version | 0.8.2 |
| Downloads | 2.1+ million total |
| License | MIT / Apache-2.0 (dual) |
| Min Rust | 1.64+ |

### Graph Types Provided

1. **Graph** - Adjacency list backed, most common use case
2. **StableGraph** - Indexes remain valid after node/edge removal
3. **GraphMap** - HashMap-backed for quick node lookup
4. **MatrixGraph** - Adjacency matrix representation
5. **Acyclic wrapper** - Enforces DAG invariants

### Algorithms Available

| Category | Algorithms |
|----------|------------|
| **Traversal** | DFS, BFS (as iterators) |
| **Shortest Path** | Dijkstra, Bellman-Ford, A*, Floyd-Warshall |
| **Minimum Spanning** | Kruskal, Prim |
| **Connectivity** | Strongly/Weakly connected components, Tarjan |
| **Isomorphism** | Subgraph isomorphism testing |
| **Flow** | Max-flow (Ford-Fulkerson variant) |
| **Matching** | Maximum matching |
| **Topological** | Topological sort, cycle detection |

### Missing Algorithms (Not in petgraph core)

- **PageRank** - Available via `graphalgs` crate
- **Centrality measures** (betweenness, closeness, eigenvector)
- **Community detection** (Louvain, Label Propagation)
- **Graph Neural Network integration**

### Optional Features

```toml
[dependencies]
petgraph = { version = "0.8", features = ["serde-1", "rayon", "dot_parser"] }
```

- `serde-1` - Serialization support
- `rayon` - Parallel iterators
- `dot_parser` - Parse GraphViz DOT format
- `std` - Standard library (disable for no_std)

### Code Example: Basic Graph Operations

```rust
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::algo::{dijkstra, is_cyclic_directed};
use petgraph::dot::{Dot, Config};

// Create a directed graph
let mut graph: DiGraph<&str, f64> = DiGraph::new();

// Add nodes
let a = graph.add_node("A");
let b = graph.add_node("B");
let c = graph.add_node("C");

// Add edges with weights
graph.add_edge(a, b, 1.5);
graph.add_edge(b, c, 2.0);
graph.add_edge(a, c, 4.0);

// Run Dijkstra from node A
let distances = dijkstra(&graph, a, None, |e| *e.weight());
println!("Distance A->C: {:?}", distances.get(&c)); // Some(3.5)

// Check for cycles
println!("Has cycle: {}", is_cyclic_directed(&graph)); // false

// Export to DOT format
println!("{:?}", Dot::with_config(&graph, &[Config::EdgeNoLabel]));
```

### Integration with RuVector

**Potential integration points:**

1. **Graph structure export** - Export RuVector's internal graph to petgraph for analysis
2. **Algorithm pipeline** - Use petgraph algorithms as preprocessing for GNN embeddings
3. **Visualization** - Export to DOT format for debugging/visualization
4. **Hybrid approach** - Use petgraph for classical algorithms, RuVector for vector similarity

**Challenges:**
- RuVector has its own graph layer (`ruvector-graph`) with Neo4j compatibility
- Different data models: petgraph is general-purpose, RuVector is vector-centric
- No direct API bridge currently documented

### Related Crates

| Crate | Purpose |
|-------|---------|
| [graphalgs](https://crates.io/crates/graphalgs) | Additional algorithms (PageRank, etc.) |
| [petgraph-gen](https://crates.io/crates/petgraph-gen) | Graph generators |
| [petgraph-graphml](https://crates.io/crates/petgraph-graphml) | GraphML export |
| [graph-api-petgraph](https://crates.io/crates/graph-api-petgraph) | Unified graph API |

---

## 2. graph crate - Alternative Implementation

**Alternatives to petgraph for specialized use cases**

There is no single "graph" crate that rivals petgraph. Instead, several alternatives exist for specific needs:

### typed_graph

| Attribute | Value |
|-----------|-------|
| Crates.io | [typed_graph](https://lib.rs/crates/typed_graph) |
| Focus | Functionality over performance |
| Use Case | Type-safe graph operations |

Focuses on compile-time type safety rather than raw performance. Good when you need strong type guarantees.

### scene-graph

| Attribute | Value |
|-----------|-------|
| Crates.io | [scene-graph](https://crates.io/crates/scene-graph) |
| Focus | Fast iteration |
| Trade-off | Slower node creation |

Similar to petgraph's StableGraph but optimized for iteration-heavy workloads (e.g., game development scene graphs).

**Key differences from petgraph:**
- Faster iteration over nodes
- Slower node creation
- No built-in algorithms (user must implement)

### gryf

| Attribute | Value |
|-----------|-------|
| GitHub | [pnevyk/gryf](https://github.com/pnevyk/gryf) |
| Focus | Lazy algorithm evaluation |
| Feature | Configurable edge directionality |

Offers lazy iterators for algorithms like topological sort. Results can be collected into `Result<Vec<VertexId>, Error>`.

### rust-graph (Deprecated)

The maintainer recommends using petgraph instead. Only maintained for bug fixes.

### Verdict

**petgraph remains the best general-purpose choice.** Use alternatives only when:
- You need type-heavy guarantees (typed_graph)
- Iteration performance is critical (scene-graph)
- You need lazy algorithm evaluation (gryf)

---

## 3. differential-dataflow - Timely Dataflow

**High-throughput, low-latency incremental computation**

| Attribute | Value |
|-----------|-------|
| Crates.io | [differential-dataflow](https://crates.io/crates/differential-dataflow) |
| Documentation | [docs.rs/differential-dataflow](https://docs.rs/differential-dataflow) |
| GitHub | [TimelyDataflow/differential-dataflow](https://github.com/TimelyDataflow/differential-dataflow) |
| Current Version | 0.12.x |
| License | MIT |
| Based On | Timely Dataflow |

### What Is It?

Differential dataflow is a data-parallel programming framework for processing large volumes of data and responding efficiently to changes in input collections. Built on [timely-dataflow](https://github.com/TimelyDataflow/timely-dataflow).

### When to Use It

| Use Case | Why Differential Dataflow? |
|----------|---------------------------|
| **Streaming analytics** | Processes incremental updates efficiently |
| **Real-time dashboards** | Sub-second latency on data changes |
| **Graph analytics at scale** | Handles billions of edges with incremental updates |
| **Iterative algorithms** | PageRank, connected components with live updates |
| **Change data capture** | React to database changes efficiently |

### Key Concepts

1. **Collection-oriented** - Transform collections with map, filter, join, reduce
2. **Incremental** - Only recomputes what changed
3. **Distributed** - Scales across threads/processes/machines
4. **Iterative** - Native support for iterate operator

### Code Example: Real-time Graph Reachability

```rust
use differential_dataflow::input::Input;
use differential_dataflow::operators::*;

fn main() {
    timely::execute_directly(|worker| {
        let (mut edges, probe) = worker.dataflow(|scope| {
            let (edges_input, edges) = scope.new_collection();

            // Define edges: (from, to)
            // Compute reachability from node 0
            let reachable = edges
                .map(|(a, b)| (a, b))
                .iterate(|reach| {
                    let edges = edges.enter(&reach.scope());
                    reach
                        .map(|(_src, dst)| dst)
                        .concat(&edges.map(|(src, _dst)| src))
                        .distinct()
                        .map(|node| (0, node))
                        .filter(|(src, _)| *src == 0)
                        .concat(&edges)
                        .distinct()
                });

            reachable.probe()
        });

        // Add edges incrementally
        edges.insert((0, 1));
        edges.insert((1, 2));
        edges.advance_to(1);
        worker.step();

        // Add more edges - only affected parts recompute
        edges.insert((2, 3));
        edges.advance_to(2);
        worker.step();
    });
}
```

### Performance Characteristics

| Scenario | Throughput |
|----------|------------|
| Batch 100k updates | ~500ms per batch |
| Single update | Sub-millisecond |
| Distributed (8 nodes) | Near-linear scaling |

### Integration with Gene Platform

**Potential uses:**
1. **Real-time knowledge graph updates** - Process streaming data source changes
2. **Incremental embedding updates** - Recompute only affected embeddings
3. **Live analytics** - Dashboard updates for connected data sources
4. **Graph algorithm streaming** - PageRank/centrality on live graphs

**Challenges:**
- Steep learning curve
- Different paradigm from batch processing
- Overkill for small datasets

---

## 4. pathfinding crate - Algorithm Focus

**Comprehensive pathfinding and graph algorithms**

| Attribute | Value |
|-----------|-------|
| Crates.io | [pathfinding](https://crates.io/crates/pathfinding) |
| Documentation | [docs.rs/pathfinding](https://docs.rs/pathfinding/) |
| Website | [rfc1149.net/devel/pathfinding](https://rfc1149.net/devel/pathfinding.html) |
| License | MIT / Apache-2.0 |
| Min Rust | 1.77.2 |

### Algorithms Beyond petgraph

| Category | pathfinding | petgraph |
|----------|-------------|----------|
| **A*** | Yes | Limited |
| **IDA*** | Yes | No |
| **Fringe Search** | Yes | No |
| **Bidirectional** | Yes | No |
| **Yen's k-shortest** | Yes | No |
| **Kuhn-Munkres** | Yes | No |
| **Grid utilities** | Yes | No |
| **Dijkstra** | Yes | Yes |
| **BFS/DFS** | Yes | Yes |
| **Connected components** | Yes | Yes |
| **Kruskal MST** | Yes | Yes |

### Key Features

1. **Grid type** - Built-in rectangular grid with auto-edges
2. **Matrix type** - Neighbor-aware data storage
3. **Generic algorithms** - Works with any type implementing traits
4. **Yen's algorithm** - Find k-shortest paths (unique to pathfinding)
5. **Hungarian algorithm** - Optimal assignment in bipartite graphs

### Code Example: A* with Grid

```rust
use pathfinding::prelude::*;

fn main() {
    // Create a 10x10 grid
    let mut grid = Grid::new(10, 10);

    // Add walkable cells
    grid.fill();

    // Add obstacles
    grid.remove_vertex((3, 3));
    grid.remove_vertex((3, 4));
    grid.remove_vertex((3, 5));

    // Find path using A*
    let start = (0, 0);
    let goal = (9, 9);

    let path = astar(
        &start,
        |&(x, y)| grid.neighbours((x, y)).map(|n| (n, 1)),
        |&(x, y)| (goal.0.abs_diff(x) + goal.1.abs_diff(y)) as usize,
        |&pos| pos == goal,
    );

    if let Some((path, cost)) = path {
        println!("Path found with cost {}: {:?}", cost, path);
    }
}
```

### Code Example: k-Shortest Paths (Yen's Algorithm)

```rust
use pathfinding::prelude::*;

fn main() {
    // Graph: A -> B -> C
    //        A -> C (direct)
    let successors = |n: &char| -> Vec<(char, usize)> {
        match n {
            'A' => vec![('B', 1), ('C', 5)],
            'B' => vec![('C', 1)],
            _ => vec![],
        }
    };

    // Find 3 shortest paths from A to C
    let paths = yen(
        &'A',
        successors,
        |&n| n == 'C',
        3,  // k = 3 paths
    );

    for (i, (path, cost)) in paths.iter().enumerate() {
        println!("Path {}: {:?} (cost: {})", i + 1, path, cost);
    }
    // Path 1: ['A', 'B', 'C'] (cost: 2)
    // Path 2: ['A', 'C'] (cost: 5)
}
```

### When to Use pathfinding Over petgraph

| Scenario | Recommended |
|----------|-------------|
| General graph operations | petgraph |
| Gaming/navigation | pathfinding |
| Grid-based problems | pathfinding |
| Need k-shortest paths | pathfinding |
| Need bipartite matching | pathfinding |
| Complex graph transformations | petgraph |
| Serialization needed | petgraph |

### Integration Strategy

Use **both** crates together:
- petgraph for graph structure and basic operations
- pathfinding for specialized algorithms when needed

```rust
// Convert petgraph to pathfinding-compatible interface
use petgraph::graph::DiGraph;
use pathfinding::prelude::astar;

fn shortest_path_in_petgraph(graph: &DiGraph<String, f64>, start: NodeIndex) {
    let path = astar(
        &start,
        |&node| graph.edges(node).map(|e| (e.target(), *e.weight() as usize)),
        |_| 0,  // No heuristic (becomes Dijkstra)
        |&node| node == goal,
    );
}
```

---

## 5. Rust Graph Databases

### 5.1 agdb - Agnesoft Graph Database

**Application-native embedded graph database with ACID compliance**

| Attribute | Value |
|-----------|-------|
| Crates.io | [agdb](https://crates.io/crates/agdb) |
| Documentation | [docs.rs/agdb](https://docs.rs/agdb/latest/agdb/) |
| GitHub | [agnesoft/agdb](https://github.com/agnesoft/agdb) |
| Current Version | 0.11.2 |
| License | Apache-2.0 |

#### Key Features

| Feature | Description |
|---------|-------------|
| **ACID Compliance** | Full ACID with write-ahead logging |
| **Memory Mapped** | Optional mmap for large datasets |
| **No Query Language** | Native object queries (no SQL/Cypher) |
| **Schema-less** | Flexible data model, no migrations |
| **Constant Time** | O(1) lookups regardless of dataset size |
| **Cluster Mode** | Enhanced consensus for consistency |

#### Storage Modes

```rust
use agdb::{Db, DbFile, DbMemory};

// Memory-mapped (default) - fastest, full ACID
let db = Db::new("./my_graph.agdb")?;

// File-based - minimal memory, full ACID, slower
let db_file = DbFile::new("./my_graph_file.agdb")?;

// In-memory - fastest, no persistence
let db_mem = DbMemory::new()?;
```

#### Code Example

```rust
use agdb::{Db, QueryBuilder};

fn main() -> Result<(), agdb::Error> {
    let mut db = Db::new("test.agdb")?;

    // Insert nodes
    db.exec_mut(&QueryBuilder::insert().nodes().aliases(&["user:1", "user:2"]).query())?;

    // Insert edge
    db.exec_mut(
        &QueryBuilder::insert()
            .edges()
            .from("user:1")
            .to("user:2")
            .values(&[("relationship", "follows").into()])
            .query()
    )?;

    // Query neighbors
    let result = db.exec(
        &QueryBuilder::select()
            .ids()
            .from(QueryBuilder::search().from("user:1").to("user:2").query())
            .query()
    )?;

    Ok(())
}
```

#### Gene Platform Fit

**Pros:**
- Constant-time lookups for data source relationships
- No ORM or query language learning curve
- Rust-native, safe, embeddable

**Cons:**
- No standard query language (Cypher, SPARQL)
- Smaller ecosystem than Neo4j/PostgreSQL
- May need custom integration layer

---

### 5.2 IndraDB - TAO-Style Graph Database

**Facebook TAO-inspired graph database**

| Attribute | Value |
|-----------|-------|
| Crates.io | [indradb](https://crates.io/crates/indradb) |
| Documentation | [docs.rs/indradb](https://docs.rs/crate/indradb/latest) |
| GitHub | [indradb/indradb](https://github.com/indradb/indradb) |
| Current Version | 4.0.0 |
| License | Apache-2.0 |
| Protocol | gRPC |

#### Key Features

| Feature | Description |
|---------|-------------|
| **TAO-inspired** | Originally modeled on Facebook's graph store |
| **Typed Edges** | Directed, typed graphs with JSON properties |
| **Multi-hop Queries** | Traverse multiple edge types |
| **Property Indexing** | Query by indexed properties |
| **Cross-Language** | gRPC bindings (Python, Rust) |
| **Pluggable Storage** | In-memory, PostgreSQL, Sled |

#### Storage Backends

| Backend | Durability | Performance | Use Case |
|---------|------------|-------------|----------|
| **Memory** | Explicit save | Fastest | Development, caching |
| **PostgreSQL** | Full | Good | Production, existing infra |
| **Sled** | Full | Fast | Embedded production |

#### Code Example

```rust
use indradb::{MemoryDatastore, Type, Vertex, SpecificVertexQuery, RangeVertexQuery};
use uuid::Uuid;

fn main() -> Result<(), indradb::Error> {
    let datastore = MemoryDatastore::default();
    let trans = datastore.transaction()?;

    // Create vertices
    let user_type = Type::new("user")?;
    let user1 = Vertex::new(user_type.clone());
    let user2 = Vertex::new(user_type);

    trans.create_vertex(&user1)?;
    trans.create_vertex(&user2)?;

    // Create edge
    let follows = Type::new("follows")?;
    trans.create_edge(&Edge::new(user1.id, follows, user2.id))?;

    // Query outbound edges
    let query = SpecificVertexQuery::single(user1.id)
        .outbound()
        .t(Type::new("follows")?);

    let edges = trans.get_edges(query)?;
    Ok(())
}
```

#### Gene Platform Fit

**Pros:**
- gRPC enables polyglot microservices
- Proven TAO model for social-graph workloads
- Multiple storage backends

**Cons:**
- Client-server model adds latency
- Less feature-rich than Neo4j
- Property queries less powerful than Cypher

---

### 5.3 Oxigraph - RDF/SPARQL Store

**Standards-compliant RDF triple store with SPARQL**

| Attribute | Value |
|-----------|-------|
| Crates.io | [oxigraph](https://crates.io/crates/oxigraph) |
| Documentation | [docs.rs/oxigraph](https://docs.rs/oxigraph) |
| GitHub | [oxigraph/oxigraph](https://github.com/oxigraph/oxigraph) |
| License | MIT / Apache-2.0 |
| Storage | RocksDB |

#### Standards Support

| Standard | Status |
|----------|--------|
| SPARQL 1.1 Query | Full |
| SPARQL 1.1 Update | Full |
| SPARQL 1.1 Federated Query | Full |
| RDF 1.2 | Behind `rdf-12` feature |
| SPARQL 1.2 | Behind `rdf-12` feature |

#### Serialization Formats

- Turtle, TriG, N-Triples, N-Quads, RDF/XML (input/output)
- SPARQL Query Results: XML, JSON, CSV, TSV

#### Code Example

```rust
use oxigraph::store::Store;
use oxigraph::model::*;
use oxigraph::sparql::QueryResults;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create persistent store
    let store = Store::open("./oxigraph_data")?;

    // Insert RDF triple
    let ex = NamedNodeRef::new("http://example.org/")?;
    let alice = NamedNodeRef::new("http://example.org/alice")?;
    let knows = NamedNodeRef::new("http://example.org/knows")?;
    let bob = NamedNodeRef::new("http://example.org/bob")?;

    store.insert(&QuadRef::new(alice, knows, bob, GraphNameRef::DefaultGraph))?;

    // SPARQL query
    let query = "SELECT ?person WHERE { <http://example.org/alice> <http://example.org/knows> ?person }";

    if let QueryResults::Solutions(solutions) = store.query(query)? {
        for solution in solutions {
            let solution = solution?;
            println!("Alice knows: {:?}", solution.get("person"));
        }
    }

    Ok(())
}
```

#### Gene Platform Fit

**Pros:**
- SPARQL is powerful for semantic queries
- Standards-compliant for interoperability
- Good for linked data / knowledge graphs

**Cons:**
- RDF model is verbose
- SPARQL learning curve
- Not optimized for vector operations

---

### 5.4 CozoDB - Datalog Graph Database

**Embedded Datalog database with vector search**

| Attribute | Value |
|-----------|-------|
| Crates.io | [cozo](https://crates.io/crates/cozo) |
| Documentation | [docs.rs/cozo](https://docs.rs/cozo) |
| GitHub | [cozodb/cozo](https://github.com/cozodb/cozo) |
| Website | [cozodb.org](https://www.cozodb.org/) |
| License | MPL-2.0 |

#### Key Features

| Feature | Description |
|---------|-------------|
| **Datalog Queries** | Powerful recursive queries |
| **Graph Algorithms** | Built-in PageRank, shortest path |
| **Vector Search** | HNSW indices since v0.6 |
| **Full-Text Search** | Added in v0.7 |
| **MinHash-LSH** | Near-duplicate detection |
| **Embeddable** | Runs in-process like SQLite |

#### Storage Backends

| Backend | Use Case |
|---------|----------|
| **RocksDB** | Production, persistent |
| **SQLite** | Lightweight, embedded |
| **Memory** | Testing, caching |
| **WASM** | Browser/mobile |

#### Performance

- 100K QPS mixed read/write/update transactions
- 250K+ QPS read-only queries
- ~50MB peak memory for 1.6M rows (RocksDB)

#### Code Example

```rust
use cozo::DbInstance;
use miette::IntoDiagnostic;

fn main() -> miette::Result<()> {
    // Create in-memory database
    let db = DbInstance::new("mem", "", Default::default()).into_diagnostic()?;

    // Create relation (table)
    db.run_script(
        ":create person {name: String => age: Int, city: String}",
        Default::default(),
    ).into_diagnostic()?;

    // Insert data
    db.run_script(
        "?[name, age, city] <- [['Alice', 30, 'NYC'], ['Bob', 25, 'LA']]
         :put person {name => age, city}",
        Default::default(),
    ).into_diagnostic()?;

    // Datalog query with recursion
    let result = db.run_script(
        "?[name, age] := *person{name, age, city}, city = 'NYC'",
        Default::default(),
    ).into_diagnostic()?;

    println!("{:?}", result);

    // Create HNSW vector index
    db.run_script(
        "::hnsw create person:embedding_idx {
            dim: 768,
            m: 16,
            ef_construction: 200,
            fields: [embedding]
        }",
        Default::default(),
    ).into_diagnostic()?;

    Ok(())
}
```

#### Gene Platform Fit

**Best candidate for Gene platform integration**

**Pros:**
- Datalog enables complex relational + graph queries
- Built-in vector search (HNSW) aligns with AI workloads
- Embeddable, cross-platform (including WASM)
- Graph algorithms included (PageRank, etc.)
- Active development

**Cons:**
- Datalog requires learning new paradigm
- Smaller community than SQL databases
- No Cypher support

---

### 5.5 Rust Graph Database Comparison

| Database | Query Language | Vector Search | Embedded | ACID | Best For |
|----------|---------------|---------------|----------|------|----------|
| **agdb** | Native Rust | No | Yes | Yes | Application-native graphs |
| **IndraDB** | Custom | No | Yes/Server | Pluggable | TAO-style social graphs |
| **Oxigraph** | SPARQL | No | Yes | No | Semantic web/RDF |
| **CozoDB** | Datalog | Yes (HNSW) | Yes | Yes | AI/ML + graph queries |
| **RuVector** | Cypher | Yes | Yes | Yes | Vector + graph hybrid |

---

## 6. RuVector's Internal Graph Usage

### Architecture Overview

RuVector Graph (`ruvector-graph`) is a high-performance graph database layer built on the core RuVector vector database with Neo4j compatibility.

| Attribute | Value |
|-----------|-------|
| Crates.io | [ruvector-graph](https://crates.io/crates/ruvector-graph) |
| Current Version | 0.1.31 |
| License | MIT |
| Documentation Coverage | 62.05% |

### Internal Dependencies

**RuVector-graph uses petgraph internally:**

```toml
[dependencies]
petgraph = "^0.6"          # Graph algorithms
roaring = "^0.10"           # Bitset operations
rayon = "^1.10"             # Parallel processing
dashmap = "^6.1"            # Concurrent hashmaps
simsimd = { version = "^5.9", optional = true }  # SIMD vectors
```

### What petgraph Provides to RuVector

1. **Graph data structures** - Adjacency lists, stable graphs
2. **Basic algorithms** - BFS, DFS, shortest paths
3. **Topological operations** - Cycle detection, ordering

### RuVector's Extensions Beyond petgraph

| Feature | petgraph | RuVector-graph |
|---------|----------|----------------|
| Cypher queries | No | Yes |
| Property graphs | Limited | Full |
| Hyperedges | No | Yes |
| Vector embeddings | No | Yes |
| GNN integration | No | Yes |
| ACID transactions | No | Yes |
| Distributed queries | No | Yes |

### API Surface

```rust
// Core types exported by ruvector-graph
pub use ruvector_graph::{
    // Database
    GraphDB,

    // Graph primitives
    Node, Edge, Hyperedge,
    Label, RelationType, Properties,

    // Transactions
    Transaction, TransactionManager,

    // Storage
    GraphStorage,
};
```

### Accessing petgraph Through RuVector

**RuVector does NOT expose petgraph directly.** However, you can:

1. **Export to petgraph format** for analysis:

```rust
use ruvector_graph::GraphDB;
use petgraph::Graph;

fn export_to_petgraph(ruvector_db: &GraphDB) -> Graph<String, f64> {
    let mut pg = Graph::new();

    // Export nodes
    let node_indices: HashMap<NodeId, NodeIndex> = ruvector_db
        .nodes()
        .map(|node| {
            let idx = pg.add_node(node.label().to_string());
            (node.id(), idx)
        })
        .collect();

    // Export edges
    for edge in ruvector_db.edges() {
        let from = node_indices[&edge.from()];
        let to = node_indices[&edge.to()];
        pg.add_edge(from, to, edge.weight());
    }

    pg
}
```

2. **Use petgraph algorithms** on exported graph:

```rust
use petgraph::algo::dijkstra;

let pg = export_to_petgraph(&ruvector_db);
let distances = dijkstra(&pg, start_node, None, |e| *e.weight());
```

### Hybrid Architecture Modules

RuVector-graph organizes functionality into:

| Module | Purpose |
|--------|---------|
| **cypher** | Query language execution |
| **hybrid** | Vector embeddings, semantic search, GNN, RAG |
| **distributed** | Sharding, replication, RPC, gossip |
| **txn** | MVCC-based transaction isolation |
| **optimizer** | Query planning and acceleration |

### Code Example: Using RuVector-graph

```rust
use ruvector_graph::{GraphDB, Node, Edge, Properties};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create database
    let db = GraphDB::open("./gene_knowledge_graph")?;

    // Create nodes with properties
    let datasource = Node::builder()
        .label("DataSource")
        .property("name", "Snowflake")
        .property("type", "warehouse")
        .build()?;

    let schema = Node::builder()
        .label("Schema")
        .property("name", "analytics")
        .build()?;

    // Create relationship
    let edge = Edge::builder()
        .from(&datasource)
        .to(&schema)
        .rel_type("CONTAINS")
        .property("created_at", chrono::Utc::now())
        .build()?;

    // Begin transaction
    let txn = db.begin_transaction()?;
    txn.create_node(&datasource)?;
    txn.create_node(&schema)?;
    txn.create_edge(&edge)?;
    txn.commit()?;

    // Cypher query
    let results = db.query(
        "MATCH (ds:DataSource)-[:CONTAINS]->(s:Schema)
         RETURN ds.name, s.name"
    )?;

    for row in results {
        println!("{} contains {}", row.get("ds.name"), row.get("s.name"));
    }

    Ok(())
}
```

### Recommendation for Gene Platform

**Use RuVector-graph as the primary graph layer** because:

1. **Already integrated** with vector search
2. **Cypher support** familiar to Neo4j users
3. **Hyperedges** for complex multi-entity relationships
4. **GNN integration** for AI-enhanced queries
5. **ACID transactions** for data integrity
6. **Distributed** capabilities for scale

**Supplement with petgraph** for:
- Offline analysis
- Visualization (DOT export)
- Algorithms not in RuVector (isomorphism, etc.)

---

## Summary and Recommendations

### Library Selection Matrix

| Use Case | Recommended Library |
|----------|-------------------|
| **Primary graph storage** | RuVector-graph |
| **Classical algorithms** | petgraph |
| **Pathfinding/gaming** | pathfinding |
| **Streaming analytics** | differential-dataflow |
| **Semantic web/RDF** | Oxigraph |
| **AI + graph queries** | CozoDB (alternative) |

### Integration Strategy for Gene

```
┌─────────────────────────────────────────────────────────┐
│                    Gene Platform                        │
├─────────────────────────────────────────────────────────┤
│  ┌──────────────────┐  ┌──────────────────────────┐     │
│  │   RuVector-graph │  │   Analysis Layer         │     │
│  │   (Primary)      │  │                          │     │
│  │                  │  │  ┌────────────────────┐  │     │
│  │  - Cypher queries│  │  │    petgraph        │  │     │
│  │  - Vector search │  │  │    (algorithms)    │  │     │
│  │  - GNN layers    │  │  └────────────────────┘  │     │
│  │  - ACID txns     │  │                          │     │
│  │  - Hyperedges    │  │  ┌────────────────────┐  │     │
│  │                  │──▶│  │   pathfinding     │  │     │
│  │  [petgraph ^0.6] │  │  │   (specialized)    │  │     │
│  │   internal       │  │  └────────────────────┘  │     │
│  └──────────────────┘  └──────────────────────────┘     │
└─────────────────────────────────────────────────────────┘
```

### Key Takeaways

1. **petgraph is foundational** - RuVector already uses it internally
2. **No need to replace RuVector-graph** - It extends petgraph with production features
3. **CozoDB is a strong alternative** if Datalog is preferred over Cypher
4. **Use pathfinding for specialized algorithms** (Yen's k-shortest, Hungarian)
5. **Consider differential-dataflow** for streaming/real-time requirements

---

## Sources

- [petgraph on crates.io](https://crates.io/crates/petgraph)
- [petgraph documentation](https://docs.rs/petgraph/)
- [petgraph GitHub](https://github.com/petgraph/petgraph)
- [pathfinding on crates.io](https://crates.io/crates/pathfinding)
- [differential-dataflow GitHub](https://github.com/TimelyDataflow/differential-dataflow)
- [agdb documentation](https://docs.rs/agdb/latest/agdb/)
- [IndraDB GitHub](https://github.com/indradb/indradb)
- [Oxigraph GitHub](https://github.com/oxigraph/oxigraph)
- [CozoDB website](https://www.cozodb.org/)
- [CozoDB GitHub](https://github.com/cozodb/cozo)
- [RuVector GitHub](https://github.com/ruvnet/ruvector)
- [ruvector-graph documentation](https://docs.rs/ruvector-graph/latest/ruvector_graph/)
