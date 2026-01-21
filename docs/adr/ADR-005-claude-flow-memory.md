# ADR-005: Claude Flow Memory Configuration

## Status
Accepted

## Context

Claude Flow provides AI agent memory capabilities including:
- Pattern storage and retrieval
- SONA learning trajectories
- Task history and outcomes
- Cross-session memory persistence

The current memory configuration needs to be fixed to properly utilize RuVector's HNSW indexing and ensure learning persists across sessions.

## Decision

Configure Claude Flow memory with:
1. **Hybrid backend** for flexibility
2. **HNSW indexing** enabled for 150x-12,500x faster pattern search
3. **Persistent storage** for learning trajectories
4. **Organized namespaces** for different memory types

### Configuration

**Hardware-Optimized for Ryzen 9 7950X3D + 187GB RAM + NVMe RAID:**

```json
{
  "memory": {
    "backend": "hybrid",
    "storagePath": "./data/claude-flow-memory",
    "hnsw": {
      "enabled": true,
      "dimensions": 768,
      "M": 48,
      "efConstruction": 400,
      "efSearch": 200
    },
    "quantization": "none",
    "cacheSize": 50000,
    "maxEntries": 10000000,
    "persistence": {
      "autoPersist": true,
      "persistInterval": 5000,
      "walMode": true
    }
  },
  "learning": {
    "sona": { "enabled": true },
    "trajectories": { "persist": true },
    "patterns": { "autoStore": true }
  },
  "ruvector": {
    "enabled": true,
    "moe": { "enabled": true, "experts": 8 },
    "ewc": { "enabled": true, "consolidationInterval": "30m" },
    "flash": { "enabled": true, "blockSize": 128 }
  }
}
```

### Namespace Layout

| Namespace | Purpose | Examples |
|-----------|---------|----------|
| `patterns` | Learned successful patterns from past tasks | "JWT auth pattern", "React form validation" |
| `trajectories` | SONA learning trajectories | Decision sequences, outcomes |
| `tasks` | Task history and outcomes | Completed tasks with results |
| `solutions` | Bug fixes and solutions for reuse | Error resolutions, workarounds |

### HNSW Parameters (Hardware-Optimized)

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `dimensions` | 768 | Matches all-mpnet-base-v2 embedding size |
| `M` | 48 | Higher connectivity - 183GB RAM supports it |
| `efConstruction` | 400 | Better index quality - NVMe handles fast writes |
| `efSearch` | 200 | Better recall - 3D V-Cache accelerates traversal |
| `cacheSize` | 50000 | Large in-memory cache - abundant RAM |
| `quantization` | none | Full f32 precision - no compression needed |

## Consequences

### Positive
- 150x-12,500x faster pattern search with HNSW
- Learning persists across sessions
- Organized memory with clear namespaces
- Consistent embedding dimensions with gene/SNP embeddings

### Negative
- Need to audit and fix current configuration
- Migration of existing memory data may be needed

### Neutral
- Storage in `./data/claude-flow-memory` directory
- Same embedding model as biomedical data (all-mpnet-base-v2)

## Verification Commands

```bash
# Test memory is working
npx @claude-flow/cli@latest memory store --key "test" --value "hello" --namespace patterns
npx @claude-flow/cli@latest memory search --query "hello" --namespace patterns
npx @claude-flow/cli@latest memory retrieve --key "test" --namespace patterns

# Check learning stats
npx @claude-flow/cli@latest hooks metrics
npx @claude-flow/cli@latest hooks intelligence stats

# Initialize memory database
npx @claude-flow/cli@latest memory init --force --verbose
```

## Implementation Steps

1. Audit current `.claude/settings.json` hooks configuration
2. Fix memory backend configuration in `claude-flow.config.json`
3. Ensure RuVector is properly initialized for Claude Flow memory operations
4. Test memory store/search/retrieve commands work correctly
5. Verify SONA learning trajectories persist across sessions
6. Configure namespaces: `patterns`, `trajectories`, `tasks`, `solutions`

## Related Decisions
- [ADR-001](./ADR-001-three-database-architecture.md): Three-Database Architecture
- [ADR-003](./ADR-003-embedding-model-selection.md): Embedding Model Selection

## References
- [Claude Flow Documentation](https://github.com/ruvnet/claude-flow)
- [HNSW Algorithm](https://arxiv.org/abs/1603.09320)
- [RuVector HNSW Implementation](https://github.com/ruvnet/ruvector)
