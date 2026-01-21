# ADR-007: Automated Learning Configuration

## Status
Accepted

## Context

Claude Flow provides multiple learning systems that are disabled by default:
- **SONA**: Self-Optimizing Neural Architecture for runtime-adaptive learning
- **EWC++**: Elastic Weight Consolidation to prevent catastrophic forgetting
- **LoRA**: Low-Rank Adaptation for efficient fine-tuning
- **MoE**: Mixture of Experts for specialized agent routing
- **Trajectories**: Recording action sequences for reinforcement learning
- **Patterns**: Storing successful approaches for reuse
- **ReasoningBank**: RETRIEVE → JUDGE → DISTILL → CONSOLIDATE pipeline

Without enabling these, Claude Flow operates in "stateless" mode with no learning between sessions.

## Decision

Enable all learning systems with hardware-optimized parameters:

### Core Learning Systems

| System | Setting | Value | Purpose |
|--------|---------|-------|---------|
| SONA | `sona.enabled` | true | Sub-millisecond adaptive learning |
| EWC++ | `ewc.enabled` | true | Prevent catastrophic forgetting |
| LoRA | `lora.enabled` | true | 128x memory compression |
| MoE | `moe.enabled` | true | Specialized expert routing |
| Trajectories | `trajectory.enabled` | true | Record action sequences |
| Patterns | `patterns.enabled` | true | Store successful approaches |
| ReasoningBank | `reasoningBank.enabled` | true | Learning pipeline |

### EWC++ Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `ewc.lambda` | 0.4 | 0.5 | Favor learning new patterns (plasticity) |
| `ewc.gamma` | 0.95 | 0.9 | Preserve important weights longer |
| `ewc.fisherSamples` | 100 | 50 | Better accuracy (CPU can handle it) |

### LoRA Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `lora.rank` | 16 | 8 | More expressive (~2MB/layer, negligible) |
| `lora.alpha` | 32 | 16 | Standard: alpha = 2 × rank |
| `lora.dropout` | 0.05 | 0.0 | Light regularization |

### MoE Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `moe.adaptiveRouting` | true | false | Learn from routing outcomes |
| `moe.experts` | 8 | 8 | 8 expert types (coder, tester, etc.) |
| `moe.topK` | 2 | 1 | Consult 2 experts per task |

### Trajectory Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `trajectory.autoRecord` | true | false | Learn from every interaction |
| `trajectory.bufferSize` | 10000 | 1000 | 10K trajectories (~100MB) |
| `trajectory.maxSteps` | 50 | 20 | Capture full task context |

### Pattern Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `patterns.autoStore` | true | false | Auto-store successful patterns |
| `patterns.clusteringEnabled` | true | false | K-means++ for pattern discovery |

### Learning Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `learning.autoLearn` | true | false | Enable automatic learning |
| `learning.learningRate` | 0.01 | 0.001 | 10x for faster convergence |
| `learning.batchSize` | 32 | 16 | Standard batch size |
| `learning.patternThreshold` | 0.85 | 0.7 | Only store high-confidence patterns |

### Hooks Parameters

| Parameter | Value | Default | Rationale |
|-----------|-------|---------|-----------|
| `hooks.trainOnSuccess` | true | false | Learn from successful operations |
| `hooks.learnFromFailures` | true | false | Learn what NOT to do |

### Background Workers

Enable 7 learning-related workers:

| Worker | Purpose |
|--------|---------|
| `ultralearn` | Deep knowledge acquisition |
| `predict` | Predictive preloading |
| `document` | Auto-documentation |
| `deepdive` | Deep code analysis |
| `refactor` | Refactoring suggestions |
| `benchmark` | Performance benchmarking |
| `preload` | Resource preloading |

## Consequences

### Positive
- Continuous improvement across sessions
- Pattern reuse reduces repeated work
- Adaptive routing to best-fit agents
- Prevents forgetting old patterns (EWC++)
- Sub-millisecond learning overhead (SONA)

### Negative
- Higher memory usage (~100MB for trajectories)
- Background workers consume CPU cycles
- More complex debugging when learning affects behavior

### Neutral
- 62 total configuration settings
- All settings in `.claude-flow/config.json`

## Configuration

```javascript
// Enable all learning
mcp__claude-flow__config_set({ key: "sona.enabled", value: true })
mcp__claude-flow__config_set({ key: "ewc.enabled", value: true })
mcp__claude-flow__config_set({ key: "lora.enabled", value: true })
mcp__claude-flow__config_set({ key: "moe.enabled", value: true })
mcp__claude-flow__config_set({ key: "trajectory.enabled", value: true })
mcp__claude-flow__config_set({ key: "patterns.enabled", value: true })
mcp__claude-flow__config_set({ key: "reasoningBank.enabled", value: true })
mcp__claude-flow__config_set({ key: "learning.autoLearn", value: true })
mcp__claude-flow__config_set({ key: "hooks.trainOnSuccess", value: true })
mcp__claude-flow__config_set({ key: "hooks.learnFromFailures", value: true })

// Apply
mcp__claude-flow__system_reset({ component: "all", confirm: true })
```

## Related Decisions
- [ADR-005](./ADR-005-claude-flow-memory.md): Claude Flow Memory Configuration
- [ADR-006](./ADR-006-hardware-optimized-hnsw.md): Hardware-Optimized HNSW Parameters

## References
- [AUTOMATED-LEARNING-GUIDE.md](../architecture/AUTOMATED-LEARNING-GUIDE.md)
- [EWC Paper](https://arxiv.org/abs/1612.00796)
- [LoRA Paper](https://arxiv.org/abs/2106.09685)
- [Claude Flow Self-Improving AI](https://github.com/ruvnet/claude-flow/issues/419)
