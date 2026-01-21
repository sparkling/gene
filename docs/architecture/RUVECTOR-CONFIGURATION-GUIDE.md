# RuVector & Claude Flow Configuration Guide

> Hardware-optimized configuration for AMD Ryzen 9 7950X3D (32 threads) with 187GB RAM

**Last Updated:** 2026-01-21
**Configuration Version:** 3.0.0-alpha
**Target Hardware:** High-performance dedicated server

## Table of Contents

1. [Hardware Specifications](#hardware-specifications)
2. [Configuration Overview](#configuration-overview)
3. [HNSW Vector Index Settings](#hnsw-vector-index-settings)
4. [Swarm & Agent Settings](#swarm--agent-settings)
5. [Memory Configuration](#memory-configuration)
6. [Neural Network Settings](#neural-network-settings)
7. [Daemon & Worker Settings](#daemon--worker-settings)
8. [Performance Settings](#performance-settings)
9. [Environment Variables](#environment-variables)
10. [Tuning Guidelines](#tuning-guidelines)

---

## Hardware Specifications

| Component | Specification |
|-----------|---------------|
| **CPU** | AMD Ryzen 9 7950X3D |
| **Cores/Threads** | 16 cores / 32 threads |
| **Max Clock** | 5.76 GHz |
| **L3 Cache** | 128 MB (3D V-Cache) |
| **RAM** | 187 GB available (196 GB total) |
| **Architecture** | x86_64, AVX-512 support |
| **Virtualization** | AMD-V enabled |

### Key Hardware Features Leveraged

- **AVX-512**: SIMD vectorization for embedding operations
- **32 threads**: Parallel agent execution and background workers
- **187GB RAM**: Large HNSW indexes, embedding caches, and memory pools
- **3D V-Cache**: Improved cache hit rates for vector search

---

## Configuration Overview

### Configuration Files

| File | Purpose |
|------|---------|
| `.mcp.json` | MCP server environment variables |
| `.claude-flow/config.json` | Runtime configuration store |
| `.claude-flow/embeddings.json` | Embedding model settings |
| `.claude-flow/daemon-state.json` | Background worker state |

### Quick Reference

```json
{
  "swarm.maxAgents": 100,
  "hnsw.m": 48,
  "hnsw.efConstruction": 400,
  "hnsw.efSearch": 200,
  "hnsw.maxElements": 100000000,
  "memory.maxEntries": 50000000,
  "neural.flashAttention": true,
  "daemon.maxConcurrent": 16,
  "performance.workerThreads": 16,
  "performance.memoryLimit": 68719476736
}
```

---

## HNSW Vector Index Settings

HNSW (Hierarchical Navigable Small Worlds) is the core algorithm for fast approximate nearest neighbor search in RuVector.

### Parameters

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `hnsw.m` | **48** | 16 | **3x default.** Higher M increases graph connectivity, improving recall. With 187GB RAM, memory overhead is acceptable. Each node maintains 48 bidirectional links. |
| `hnsw.efConstruction` | **400** | 100 | **4x default.** Controls build-time accuracy. Higher values find more optimal connections during index creation. Slower build time but better index quality. |
| `hnsw.efSearch` | **200** | 50 | **4x default.** Size of candidate queue during search. Higher values improve recall at cost of latency. With fast CPU, latency impact is minimal. |
| `hnsw.maxElements` | **100,000,000** | 10,000,000 | **10x default.** Maximum vectors the index can hold. With 187GB RAM, we can support 100M vectors (384-dim embeddings ≈ 150GB at peak). |

### Trade-off Analysis

```
┌─────────────────────────────────────────────────────────────────┐
│                    HNSW Parameter Trade-offs                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  M (Connections)     ────────────────────────────────────────►  │
│  Low (8)                    Medium (16-32)              High (48)│
│  • Fast build               • Balanced                • Slow build│
│  • Low memory               • Good recall             • High memory│
│  • Poor recall              • Moderate speed          • Best recall│
│                                                                  │
│  efConstruction      ────────────────────────────────────────►  │
│  Low (100)                  Medium (200)             High (400) │
│  • Fast indexing            • Good quality           • Slow index │
│  • Lower quality            • Balanced               • Best quality│
│  • Quick updates            • Normal updates         • Rare updates│
│                                                                  │
│  efSearch            ────────────────────────────────────────►  │
│  Low (50)                   Medium (100)             High (200) │
│  • Fast queries             • Balanced               • Slow queries│
│  • Lower recall             • Good recall            • Best recall │
│  • Real-time apps           • General use            • Accuracy-first│
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### Why These Values?

1. **m=48**: With 187GB RAM, we prioritize recall over memory. Each additional connection costs ~8 bytes per vector. For 100M vectors: 100M × 48 × 8 = 38.4GB, well within budget.

2. **efConstruction=400**: Index builds are infrequent (data ingestion). We accept slower builds for permanently better index quality.

3. **efSearch=200**: Search is the hot path. Our 7950X3D's cache and clock speed minimize latency impact from larger candidate sets.

### Performance Expectations

| Vector Count | Search Latency | Memory Usage |
|--------------|----------------|--------------|
| 1,000 | <0.2ms | ~2MB |
| 100,000 | <0.5ms | ~200MB |
| 1,000,000 | <1ms | ~2GB |
| 10,000,000 | <2ms | ~20GB |
| 100,000,000 | <5ms | ~150GB |

---

## Swarm & Agent Settings

### Parameters

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `swarm.maxAgents` | **100** | 8 | **12.5x default.** With 32 threads and 187GB RAM, we can run many concurrent agents. Each agent uses ~50-200MB RAM. |
| `swarm.autoScale.minAgents` | **4** | 1 | Always maintain 4 ready agents for instant task pickup. Reduces cold-start latency. |
| `swarm.autoScale.maxAgents` | **100** | 8 | Match maxAgents for full auto-scaling capability. |
| `swarm.topology` | **hierarchical-mesh** | hierarchical | Combines queen coordination with peer-to-peer communication. Best for 10+ agents. |
| `swarm.autoScale` | **true** | true | Enable dynamic scaling based on workload. |

### Topology Selection

```
┌─────────────────────────────────────────────────────────────────┐
│                     Topology Comparison                          │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  hierarchical          mesh              hierarchical-mesh       │
│       ┌─┐              ┌─┐───┌─┐              ┌─┐                │
│       │Q│              │ │   │ │              │Q│                │
│       └┬┘              └─┼───┼─┘              └┬┘                │
│     ┌──┼──┐            ┌─┼───┼─┐           ┌──┼──┐               │
│     │  │  │            │ │   │ │           │  │  │               │
│    ┌┴┐┌┴┐┌┴┐          ┌┴─┴┐ ┌┴─┴┐        ┌─┴──┴──┴─┐             │
│    │W││W││W│          │   │ │   │        │ W───W───W│             │
│    └─┘└─┘└─┘          └───┘ └───┘        └─────────┘             │
│                                                                  │
│  • Queen control      • Peer-to-peer    • Queen + peer comm     │
│  • Anti-drift         • Resilient       • Best scalability      │
│  • Best for <10       • Best for chaos  • Best for 10-100       │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### Why hierarchical-mesh?

- **Queen coordination**: Prevents agent drift and conflicting decisions
- **Mesh communication**: Workers can share context without routing through queen
- **Scalability**: Handles 100 agents efficiently
- **Fault tolerance**: Mesh links provide redundancy

---

## Memory Configuration

### Parameters

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `memory.maxEntries` | **50,000,000** | 1,000,000 | **50x default.** With 187GB RAM, we can store 50M memory entries. Average entry ~500 bytes = 25GB. |
| `memory.persistInterval` | **30,000** (30s) | 60,000 | **2x faster.** More frequent persistence reduces data loss risk. Fast NVMe makes this cheap. |
| `embeddings.cacheSize` | **8,192** | 256 | **32x default.** Cache 8K embeddings in memory. Each 384-dim embedding = 1.5KB. Total: 12MB cache. |

### Memory Budget Allocation

```
┌─────────────────────────────────────────────────────────────────┐
│              Memory Budget (187GB Total)                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌─────────────────────────────────────────────────────────┐    │
│  │ HNSW Index (peak)                              150 GB   │    │
│  └─────────────────────────────────────────────────────────┘    │
│  ┌───────────────────────────┐                                  │
│  │ Memory Store               25 GB                        │    │
│  └───────────────────────────┘                                  │
│  ┌─────────────┐                                                │
│  │ Agents       8 GB (100 × 80MB avg)                      │    │
│  └─────────────┘                                                │
│  ┌──────┐                                                       │
│  │ Cache 0.5 GB                                            │    │
│  └──────┘                                                       │
│  ┌───┐                                                          │
│  │OS  4 GB                                                 │    │
│  └───┘                                                          │
│                                                                  │
│  Reserved headroom: ~0 GB (peak usage scenario)                 │
│  Typical usage: 40-60 GB                                        │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Neural Network Settings

### Parameters

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `neural.enabled` | **true** | true | Enable neural pattern learning and routing. |
| `neural.flashAttention` | **true** | false | **Enable 2.49-7.47x speedup.** Flash Attention reduces memory from O(N²) to O(N) and improves throughput. |
| `neural.epochs` | **100** | 50 | **2x default.** More training epochs for better pattern convergence. Fast CPU makes this feasible. |
| `lora.rank` | **16** | 8 | **2x default.** Higher LoRA rank captures more nuanced adaptations. Memory cost: negligible. |

### Flash Attention Benefits

```
┌─────────────────────────────────────────────────────────────────┐
│                  Flash Attention vs Standard                     │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Memory Usage (sequence length 8192):                           │
│  ├── Standard: O(N²) = 256 MB                                   │
│  └── Flash:    O(N)  = 32 MB   (8x reduction)                   │
│                                                                  │
│  Throughput (tokens/second):                                    │
│  ├── Standard: 1,000 tokens/s                                   │
│  └── Flash:    2,490-7,470 tokens/s  (2.49-7.47x faster)        │
│                                                                  │
│  Why enable?                                                    │
│  • Reduces GPU/CPU memory pressure                              │
│  • Enables longer context windows                               │
│  • IO-aware computation maximizes cache utilization             │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

### LoRA Rank Selection

| Rank | Adaptation Quality | Memory Overhead |
|------|-------------------|-----------------|
| 4 | Basic patterns | Minimal |
| 8 | Good patterns (default) | Low |
| **16** | **Fine-grained patterns** | **Moderate** |
| 32 | Maximum adaptation | Higher |

We chose rank=16 because:
- 2x capacity for pattern capture
- Memory overhead is ~2MB per adapted layer
- Better specialization for domain-specific tasks

---

## Daemon & Worker Settings

### Parameters

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `daemon.maxConcurrent` | **16** | 2 | **8x default.** Half of 32 threads dedicated to background workers. Leaves 16 threads for main workload. |
| `daemon.resourceThresholds.maxCpuLoad` | **24** | 2 | **12x default.** Allow workers to use up to 75% of CPU (24 of 32 threads) before throttling. |

### Background Worker Types

| Worker | Interval | Priority | Description |
|--------|----------|----------|-------------|
| `map` | 15 min | Normal | Codebase structure mapping |
| `audit` | 10 min | Critical | Security vulnerability scanning |
| `optimize` | 15 min | High | Performance optimization suggestions |
| `consolidate` | 30 min | Low | Memory consolidation (EWC++) |
| `testgaps` | 20 min | Normal | Test coverage analysis |
| `predict` | 10 min | Low | Predictive preloading |
| `document` | 60 min | Low | Auto-documentation generation |

### Worker Thread Allocation

```
┌─────────────────────────────────────────────────────────────────┐
│              Thread Allocation (32 Total)                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Main Workload (Agents, Tasks)                                  │
│  ████████████████░░░░░░░░░░░░░░░░  16 threads (50%)             │
│                                                                  │
│  Background Workers (Daemon)                                    │
│  ░░░░░░░░░░░░░░░░████████████████  16 threads (50%)             │
│                                                                  │
│  CPU Load Threshold: 24/32 (75%)                                │
│  ████████████████████████░░░░░░░░  24 threads max               │
│                                                                  │
│  Reserved for OS/System                                         │
│  ░░░░░░░░░░░░░░░░░░░░░░░░████████   8 threads (25%)             │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Performance Settings

### Parameters

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `performance.workerThreads` | **16** | 4 | **4x default.** Half of available threads for Node.js worker pool. |
| `performance.memoryLimit` | **68,719,476,736** (64GB) | 4GB | **16x default.** Allocate 64GB to Claude Flow, leaving 123GB for OS and other processes. |

### Node.js Heap Configuration

```bash
NODE_OPTIONS="--max-old-space-size=65536"  # 64GB heap
```

**Rationale:**
- Default Node.js heap: 4GB
- Our allocation: 64GB (34% of 187GB)
- Prevents OOM errors during large vector operations
- Leaves 123GB for HNSW mmap and OS cache

---

## Environment Variables

Complete `.mcp.json` environment configuration:

```json
{
  "mcpServers": {
    "claude-flow": {
      "command": "npx",
      "args": ["claude-flow@v3alpha", "mcp", "start"],
      "env": {
        "CLAUDE_FLOW_MODE": "v3",
        "CLAUDE_FLOW_HOOKS_ENABLED": "true",
        "CLAUDE_FLOW_TOPOLOGY": "hierarchical-mesh",
        "CLAUDE_FLOW_MAX_AGENTS": "100",
        "CLAUDE_FLOW_MEMORY_BACKEND": "hybrid",
        "CLAUDE_FLOW_MEMORY_PATH": ".swarm/memory.db",
        "CLAUDE_FLOW_HNSW_ENABLED": "true",
        "CLAUDE_FLOW_HNSW_M": "48",
        "CLAUDE_FLOW_HNSW_EF_CONSTRUCTION": "400",
        "CLAUDE_FLOW_HNSW_EF_SEARCH": "200",
        "CLAUDE_FLOW_FLASH_ATTENTION": "true",
        "CLAUDE_FLOW_WORKER_THREADS": "16",
        "CLAUDE_FLOW_MEMORY_MAX_ENTRIES": "50000000",
        "CLAUDE_FLOW_EMBEDDINGS_CACHE_SIZE": "8192",
        "CLAUDE_FLOW_DAEMON_MAX_CONCURRENT": "16",
        "CLAUDE_FLOW_NEURAL_ENABLED": "true",
        "CLAUDE_FLOW_LORA_RANK": "16",
        "NODE_OPTIONS": "--max-old-space-size=65536"
      }
    }
  }
}
```

---

## Tuning Guidelines

### For Different Hardware Profiles

#### Mid-Range Server (8 cores, 32GB RAM)

```json
{
  "swarm.maxAgents": 20,
  "hnsw.m": 24,
  "hnsw.efConstruction": 200,
  "hnsw.efSearch": 100,
  "hnsw.maxElements": 10000000,
  "memory.maxEntries": 5000000,
  "daemon.maxConcurrent": 4,
  "performance.workerThreads": 4,
  "NODE_OPTIONS": "--max-old-space-size=16384"
}
```

#### Entry-Level (4 cores, 16GB RAM)

```json
{
  "swarm.maxAgents": 8,
  "hnsw.m": 16,
  "hnsw.efConstruction": 100,
  "hnsw.efSearch": 50,
  "hnsw.maxElements": 1000000,
  "memory.maxEntries": 1000000,
  "daemon.maxConcurrent": 2,
  "performance.workerThreads": 2,
  "NODE_OPTIONS": "--max-old-space-size=8192"
}
```

### Scaling Formulas

```
maxAgents = min(100, threads × 3)
hnsw.m = min(64, 16 + (RAM_GB / 8))
hnsw.efConstruction = min(500, 100 + (RAM_GB / 2))
hnsw.efSearch = min(300, 50 + (RAM_GB / 4))
daemon.maxConcurrent = max(2, threads / 2)
workerThreads = max(2, threads / 2)
memoryLimit = RAM_GB × 0.34 × 1024³
```

### Monitoring & Adjustment

```bash
# Check current performance
npx @claude-flow/cli@latest performance metrics --metric all

# Run benchmarks
npx @claude-flow/cli@latest performance benchmark --suite all

# Detect bottlenecks
npx @claude-flow/cli@latest performance bottleneck --deep true
```

---

## References

- [OpenSearch HNSW Guide](https://opensearch.org/blog/a-practical-guide-to-selecting-hnsw-hyperparameters/)
- [Pinecone HNSW Explanation](https://www.pinecone.io/learn/series/faiss/hnsw/)
- [RuVector npm Documentation](https://www.npmjs.com/package/ruvector)
- [Claude Flow GitHub](https://github.com/ruvnet/claude-flow)
- [Flash Attention Paper](https://arxiv.org/abs/2205.14135)

---

*Configuration applied: 2026-01-21*
*Hardware: AMD Ryzen 9 7950X3D, 187GB RAM*
