# Claude Flow Configuration Summary (ADR 005-008)

This document summarizes the Claude Flow configuration decisions for the Gene platform.

## Overview

Claude Flow is configured for high-performance operation on the Gene platform hardware:
- **CPU**: AMD Ryzen 9 7950X3D (16C/32T, 128MB 3D V-Cache)
- **RAM**: 187GB total
- **Storage**: 2x 1.7TB NVMe RAID

---

## ADR-005: Claude Flow Memory Configuration

**Status**: Accepted

**Decision**: Configure Claude Flow memory with hybrid backend, HNSW indexing, and organized namespaces.

**Configuration Highlights**:
```json
{
  "memory": {
    "backend": "hybrid",
    "hnsw": { "enabled": true, "dimensions": 768, "M": 48 },
    "cacheSize": 50000,
    "persistence": { "autoPersist": true }
  },
  "learning": {
    "sona": { "enabled": true },
    "trajectories": { "persist": true },
    "patterns": { "autoStore": true }
  }
}
```

**Namespace Layout**:
| Namespace | Purpose |
|-----------|---------|
| `patterns` | Learned successful patterns from past tasks |
| `trajectories` | SONA learning trajectories |
| `tasks` | Task history and outcomes |
| `solutions` | Bug fixes and solutions for reuse |

**Benefits**:
- 150x-12,500x faster pattern search with HNSW
- Learning persists across sessions
- Consistent 768-dimension embeddings

---

## ADR-006: Hardware-Optimized HNSW Parameters

**Status**: Accepted

**Decision**: Adopt aggressive HNSW parameters optimized for the available hardware.

**HNSW Settings**:
| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `M` | 48 | 16 | More connections = better recall; 183GB RAM abundant |
| `efConstruction` | 400 | 200 | Higher quality index; NVMe handles fast writes |
| `efSearch` | 200 | 100 | Better recall; 3D V-Cache accelerates traversal |
| `maxElements` | 100M | 1M | NVMe + RAM supports 100M vectors |

**Memory Settings**:
| Parameter | Value | Default |
|-----------|-------|---------|
| `maxEntries` | 50M | 10K |
| `embeddings.cacheSize` | 8,192 | 256 |
| `performance.memoryLimit` | 64GB | 4GB |
| `performance.workerThreads` | 16 | 4 |

**Performance**:
- 150x-12,500x faster similarity search vs brute force
- Better recall (99%+ at k=10)
- No quantization needed (full f32 precision)

---

## ADR-007: Automated Learning Configuration

**Status**: Accepted

**Decision**: Enable all learning systems with hardware-optimized parameters.

**Learning Systems Enabled**:
| System | Purpose |
|--------|---------|
| SONA | Sub-millisecond adaptive learning |
| EWC++ | Prevent catastrophic forgetting |
| LoRA | 128x memory compression |
| MoE | Specialized expert routing (8 experts) |
| Trajectories | Record action sequences |
| Patterns | Store successful approaches |
| ReasoningBank | RETRIEVE - JUDGE - DISTILL - CONSOLIDATE pipeline |

**Key Parameters**:
| System | Parameter | Value |
|--------|-----------|-------|
| EWC++ | lambda | 0.4 (favor learning) |
| EWC++ | gamma | 0.95 (preserve weights) |
| LoRA | rank | 16 |
| MoE | topK | 2 (consult 2 experts) |
| Trajectory | bufferSize | 10,000 |
| Learning | learningRate | 0.01 (10x default) |

**Background Workers Enabled**:
- ultralearn, predict, document, deepdive, refactor, benchmark, preload

---

## ADR-008: MCP Server Management

**Status**: Accepted

**Decision**: Use `system_reset` as the primary restart method; `config_set` for runtime changes.

**Restart Methods**:
| Method | Command | Reloads Config | Preserves Session |
|--------|---------|----------------|-------------------|
| System Reset | `system_reset({ component: "all", confirm: true })` | Yes | Yes |
| MCP Reconnect | `/mcp reconnect claude-flow` | No | Yes |
| Full Restart | Exit + relaunch Claude Code | Yes | No |

**Configuration Sources (Priority Order)**:
1. MCP tool calls (immediate, session-scoped)
2. `.claude-flow/config.json` (persistent, reloaded on system_reset)
3. `.mcp.json` env vars (startup only, requires full restart)

**Common Operations**:
```javascript
// Set runtime config
mcp__claude-flow__config_set({ key: "hnsw.m", value: 48 })

// Apply changes
mcp__claude-flow__system_reset({ component: "all", confirm: true })

// Verify
mcp__claude-flow__system_status({ verbose: true })
```

---

## Quick Reference

### Verify Learning Active
```bash
npx @claude-flow/cli@latest hooks metrics
npx @claude-flow/cli@latest hooks intelligence stats
```

### Test Memory
```bash
npx @claude-flow/cli@latest memory store --key "test" --value "hello" --namespace patterns
npx @claude-flow/cli@latest memory search --query "hello"
```

### Check System Status
```javascript
mcp__claude-flow__system_status({ verbose: true })
mcp__claude-flow__hooks_intelligence({ showStatus: true })
```

---

## Related Documents

- [ADR-005](../adr/adr-005-claude-flow-memory.md): Full memory configuration decision
- [ADR-006](../adr/adr-006-hardware-optimized-hnsw.md): HNSW optimization details
- [ADR-007](../adr/adr-007-automated-learning-configuration.md): Learning configuration details
- [ADR-008](../adr/adr-008-mcp-server-management.md): MCP server management details
- [Automated Learning Guide](../architecture/automated-learning-guide.md): Implementation guide
- [MCP Configuration Guide](../architecture/mcp-configuration-guide.md): MCP setup guide
