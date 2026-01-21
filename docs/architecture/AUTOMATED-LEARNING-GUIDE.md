# Claude Flow Automated Learning Configuration Guide

> Complete reference for self-learning, neural adaptation, and pattern recognition settings

**Last Updated:** 2026-01-21
**Total Settings:** 62 configurations
**Target Hardware:** AMD Ryzen 9 7950X3D (32 threads), 187GB RAM

## Table of Contents

1. [Quick Start](#quick-start)
2. [How to Set Configuration Values](#how-to-set-configuration-values)
3. [How to Restart Claude Flow](#how-to-restart-claude-flow)
4. [Learning Pipeline Overview](#learning-pipeline-overview)
5. [SONA Settings](#sona-settings)
6. [EWC++ Settings](#ewc-settings)
7. [LoRA Settings](#lora-settings)
8. [MoE Settings](#moe-settings)
9. [Trajectory Settings](#trajectory-settings)
10. [Pattern Settings](#pattern-settings)
11. [ReasoningBank Settings](#reasoningbank-settings)
12. [Learning Settings](#learning-settings)
13. [Hooks Settings](#hooks-settings)
14. [Background Workers](#background-workers)
15. [Complete Configuration Reference](#complete-configuration-reference)

---

## Quick Start

### Enable All Learning Features

```bash
# Via MCP tools (recommended)
mcp__claude-flow__config_set({ key: "learning.autoLearn", value: true })
mcp__claude-flow__config_set({ key: "sona.enabled", value: true })
mcp__claude-flow__config_set({ key: "ewc.enabled", value: true })
mcp__claude-flow__config_set({ key: "lora.enabled", value: true })

# Restart to apply
mcp__claude-flow__system_reset({ component: "all", confirm: true })
```

### Verify Learning Status

```bash
# Check intelligence system
mcp__claude-flow__hooks_intelligence({ showStatus: true })

# Check learning statistics
mcp__claude-flow__hooks_intelligence_stats({ detailed: true })
```

---

## How to Set Configuration Values

### Method 1: MCP Tools (Recommended)

Use the `config_set` MCP tool for immediate application:

```javascript
// Set a single value
mcp__claude-flow__config_set({
  key: "learning.autoLearn",
  value: true
})

// Set numeric values
mcp__claude-flow__config_set({
  key: "ewc.lambda",
  value: 0.4
})

// Set nested values
mcp__claude-flow__config_set({
  key: "daemon.workers.ultralearn.enabled",
  value: true
})
```

### Method 2: CLI Commands

```bash
# Set via CLI
npx @claude-flow/cli@latest config set learning.autoLearn true
npx @claude-flow/cli@latest config set ewc.lambda 0.4

# List all settings
npx @claude-flow/cli@latest config list

# Get specific setting
npx @claude-flow/cli@latest config get learning.autoLearn
```

### Method 3: Direct File Edit

Edit `.claude-flow/config.json`:

```json
{
  "values": {
    "learning.autoLearn": true,
    "sona.enabled": true,
    "ewc.enabled": true,
    "ewc.lambda": 0.4
  }
}
```

### Method 4: Environment Variables

Edit `.mcp.json` for MCP server settings:

```json
{
  "mcpServers": {
    "claude-flow": {
      "env": {
        "CLAUDE_FLOW_LEARNING_AUTO": "true",
        "CLAUDE_FLOW_SONA_ENABLED": "true",
        "CLAUDE_FLOW_EWC_ENABLED": "true"
      }
    }
  }
}
```

---

## How to Restart Claude Flow

### Method 1: System Reset (Recommended - No Session Loss)

Resets internal state without losing MCP connection:

```javascript
// Via MCP tool
mcp__claude-flow__system_reset({
  component: "all",  // Options: "all", "metrics", "agents", "tasks"
  confirm: true
})
```

```bash
# Via CLI
npx @claude-flow/cli@latest system reset --component all --confirm
```

**Result:** Uptime resets to 0, all configs reloaded, no session loss.

### Method 2: Daemon Restart

Restart background workers:

```bash
# Stop daemon
npx @claude-flow/cli@latest daemon stop

# Start daemon
npx @claude-flow/cli@latest daemon start
```

### Method 3: MCP Reconnect (Session Preserved)

Use Claude Code's built-in reconnect:

```
/mcp reconnect claude-flow
```

### Method 4: Full Restart (Session Lost)

For major configuration changes (e.g., `.mcp.json` environment variables):

1. Exit Claude Code completely
2. Relaunch Claude Code
3. MCP servers reinitialize automatically

### Method 5: Swarm Restart

Reset swarm without affecting other components:

```javascript
mcp__claude-flow__swarm_shutdown({ graceful: true })
mcp__claude-flow__swarm_init({ topology: "hierarchical-mesh", maxAgents: 100 })
```

### Verification After Restart

```javascript
// Check uptime (should be near 0)
mcp__claude-flow__system_status({ verbose: true })

// Verify learning systems active
mcp__claude-flow__hooks_intelligence({ showStatus: true })
```

---

## Learning Pipeline Overview

Claude Flow implements a 4-stage learning pipeline:

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    AUTOMATED LEARNING PIPELINE                           │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌────────┐ │
│  │   RETRIEVE   │───▶│    JUDGE     │───▶│   DISTILL    │───▶│CONSOLI-│ │
│  │    (HNSW)    │    │  (Verdict)   │    │   (LoRA)     │    │ DATE   │ │
│  └──────────────┘    └──────────────┘    └──────────────┘    │(EWC++) │ │
│        │                    │                   │             └────────┘ │
│        ▼                    ▼                   ▼                  │     │
│  Fetch similar       Rate success/       Extract key         Prevent    │
│  patterns via        failure with        learnings via       forgetting │
│  vector search       verdicts            compression         of old     │
│                                                               patterns   │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────────┐    │
│  │                        SONA ORCHESTRATOR                         │    │
│  │         Sub-millisecond adaptive learning (<0.05ms)             │    │
│  └─────────────────────────────────────────────────────────────────┘    │
│                                                                          │
│  ┌─────────────────────────────────────────────────────────────────┐    │
│  │                     MoE EXPERT ROUTING                           │    │
│  │    8 experts: coder, tester, reviewer, architect, security,     │    │
│  │              performance, researcher, coordinator               │    │
│  └─────────────────────────────────────────────────────────────────┘    │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## SONA Settings

**SONA** (Self-Optimizing Neural Architecture) enables runtime-adaptive learning without expensive retraining.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `sona.enabled` | **true** | false | Enable the SONA optimizer for sub-millisecond learning |

### How SONA Works

1. Records trajectories (sequences of actions)
2. Evaluates success/failure verdicts
3. Adapts routing decisions in real-time
4. Learning overhead: <0.05ms per decision

### Set SONA Settings

```javascript
mcp__claude-flow__config_set({ key: "sona.enabled", value: true })
```

---

## EWC++ Settings

**EWC++** (Elastic Weight Consolidation) prevents catastrophic forgetting when learning new patterns.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `ewc.enabled` | **true** | false | Enable EWC++ consolidation |
| `ewc.lambda` | **0.4** | 0.5 | Constraint strength (0.0-1.0). Lower = more plasticity, higher = more stability |
| `ewc.gamma` | **0.95** | 0.9 | Fisher decay rate. Higher = longer memory of important weights |
| `ewc.fisherSamples` | **100** | 50 | Samples for Fisher information estimation. More = better accuracy |

### Why These Values?

- **lambda=0.4**: Balanced between learning new patterns (plasticity) and retaining old ones (stability). We chose slightly lower than default (0.5) to favor learning new patterns on this high-performance system.

- **gamma=0.95**: High decay rate preserves important weight information longer. Good for systems that accumulate knowledge over time.

- **fisherSamples=100**: 2x default for better estimation accuracy. Our CPU can handle the extra computation easily.

### Set EWC Settings

```javascript
mcp__claude-flow__config_set({ key: "ewc.enabled", value: true })
mcp__claude-flow__config_set({ key: "ewc.lambda", value: 0.4 })
mcp__claude-flow__config_set({ key: "ewc.gamma", value: 0.95 })
mcp__claude-flow__config_set({ key: "ewc.fisherSamples", value: 100 })
```

---

## LoRA Settings

**LoRA** (Low-Rank Adaptation) enables efficient fine-tuning with 128x memory compression.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `lora.enabled` | **true** | false | Enable LoRA adapters |
| `lora.rank` | **16** | 8 | Rank of adaptation matrices. Higher = more expressive, more memory |
| `lora.alpha` | **32** | 16 | Scaling factor. Typically 2x rank |
| `lora.dropout` | **0.05** | 0.0 | Regularization to prevent overfitting |

### Why These Values?

- **rank=16**: 2x default for more expressive adaptations. Memory cost is ~2MB per adapted layer, negligible with 187GB RAM.

- **alpha=32**: Standard practice is alpha = 2 × rank. Controls the magnitude of adaptation.

- **dropout=0.05**: Light regularization prevents overfitting to recent patterns while still learning effectively.

### Rank Selection Guide

| Rank | Use Case | Memory per Layer |
|------|----------|------------------|
| 4 | Simple patterns, low memory | ~0.5 MB |
| 8 | Default, balanced | ~1 MB |
| **16** | **Fine-grained patterns (our choice)** | **~2 MB** |
| 32 | Maximum adaptation | ~4 MB |

### Set LoRA Settings

```javascript
mcp__claude-flow__config_set({ key: "lora.enabled", value: true })
mcp__claude-flow__config_set({ key: "lora.rank", value: 16 })
mcp__claude-flow__config_set({ key: "lora.alpha", value: 32 })
mcp__claude-flow__config_set({ key: "lora.dropout", value: 0.05 })
```

---

## MoE Settings

**MoE** (Mixture of Experts) routes tasks to specialized expert agents.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `moe.enabled` | **true** | false | Enable MoE routing |
| `moe.adaptiveRouting` | **true** | false | Learn from routing outcomes |
| `moe.experts` | **8** | 8 | Number of expert types |
| `moe.topK` | **2** | 1 | Number of experts consulted per task |

### Expert Types

1. **coder** - Implementation tasks
2. **tester** - Test writing and validation
3. **reviewer** - Code review and quality
4. **architect** - System design decisions
5. **security** - Security analysis
6. **performance** - Optimization tasks
7. **researcher** - Information gathering
8. **coordinator** - Multi-agent orchestration

### Why topK=2?

- Consults 2 experts per task instead of 1
- Provides diverse perspectives
- Slight increase in computation, but better results
- Our 32-thread CPU handles parallel expert evaluation easily

### Set MoE Settings

```javascript
mcp__claude-flow__config_set({ key: "moe.enabled", value: true })
mcp__claude-flow__config_set({ key: "moe.adaptiveRouting", value: true })
mcp__claude-flow__config_set({ key: "moe.experts", value: 8 })
mcp__claude-flow__config_set({ key: "moe.topK", value: 2 })
```

---

## Trajectory Settings

**Trajectories** record sequences of actions for SONA/ReasoningBank learning.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `trajectory.enabled` | **true** | false | Enable trajectory recording |
| `trajectory.autoRecord` | **true** | false | Automatically record all sessions |
| `trajectory.bufferSize` | **10000** | 1000 | Max trajectories in memory |
| `trajectory.maxSteps` | **50** | 20 | Max steps per trajectory |

### Why These Values?

- **bufferSize=10000**: With 187GB RAM, we can store 10K trajectories easily (~100MB). More data = better learning.

- **maxSteps=50**: Most meaningful task sequences complete in <50 steps. Captures full context without excessive storage.

- **autoRecord=true**: Learns from every interaction without manual intervention.

### Trajectory Structure

```json
{
  "trajectoryId": "traj_123",
  "goal": "Implement authentication",
  "steps": [
    { "action": "read", "input": "auth.ts", "result": "success" },
    { "action": "edit", "input": "add JWT handler", "result": "success" },
    { "action": "test", "input": "run tests", "result": "pass" }
  ],
  "verdict": "success",
  "timestamp": "2026-01-21T23:00:00Z"
}
```

### Set Trajectory Settings

```javascript
mcp__claude-flow__config_set({ key: "trajectory.enabled", value: true })
mcp__claude-flow__config_set({ key: "trajectory.autoRecord", value: true })
mcp__claude-flow__config_set({ key: "trajectory.bufferSize", value: 10000 })
mcp__claude-flow__config_set({ key: "trajectory.maxSteps", value: 50 })
```

---

## Pattern Settings

**Patterns** store successful approaches for reuse via HNSW vector search.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `patterns.enabled` | **true** | false | Enable pattern storage |
| `patterns.autoStore` | **true** | false | Auto-store successful patterns |
| `patterns.clusteringEnabled` | **true** | false | K-means++ clustering for pattern discovery |

### Why Enable Clustering?

K-means++ clustering groups similar patterns, enabling:
- Faster pattern retrieval (search cluster centroids first)
- Discovery of common approaches across tasks
- Anomaly detection (patterns far from clusters)

### Set Pattern Settings

```javascript
mcp__claude-flow__config_set({ key: "patterns.enabled", value: true })
mcp__claude-flow__config_set({ key: "patterns.autoStore", value: true })
mcp__claude-flow__config_set({ key: "patterns.clusteringEnabled", value: true })
```

---

## ReasoningBank Settings

**ReasoningBank** stores reasoning traces and verdict judgments.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `reasoningBank.enabled` | **true** | false | Enable ReasoningBank |
| `reasoningBank.verdictTracking` | **true** | false | Track success/failure verdicts |

### The RETRIEVE → JUDGE → DISTILL → CONSOLIDATE Pipeline

1. **RETRIEVE**: HNSW searches for similar past situations (150x-12,500x faster)
2. **JUDGE**: Verdict system rates success/failure
3. **DISTILL**: LoRA extracts key learnings (128x compression)
4. **CONSOLIDATE**: EWC++ prevents forgetting old patterns

### Set ReasoningBank Settings

```javascript
mcp__claude-flow__config_set({ key: "reasoningBank.enabled", value: true })
mcp__claude-flow__config_set({ key: "reasoningBank.verdictTracking", value: true })
```

---

## Learning Settings

Core learning hyperparameters.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `learning.autoLearn` | **true** | false | Enable automatic learning |
| `learning.learningRate` | **0.01** | 0.001 | Gradient descent speed |
| `learning.batchSize` | **32** | 16 | Samples per training iteration |
| `learning.patternThreshold` | **0.85** | 0.7 | Minimum confidence to store pattern |

### Why These Values?

- **learningRate=0.01**: 10x default for faster learning. Our system sees enough examples to converge quickly without instability.

- **batchSize=32**: Standard batch size for neural training. Balances gradient estimation quality with memory usage.

- **patternThreshold=0.85**: Only store high-confidence patterns. Prevents noise from low-quality examples.

### Set Learning Settings

```javascript
mcp__claude-flow__config_set({ key: "learning.autoLearn", value: true })
mcp__claude-flow__config_set({ key: "learning.learningRate", value: 0.01 })
mcp__claude-flow__config_set({ key: "learning.batchSize", value: 32 })
mcp__claude-flow__config_set({ key: "learning.patternThreshold", value: 0.85 })
```

---

## Hooks Settings

Hooks trigger learning on specific events.

| Setting | Value | Default | Rationale |
|---------|-------|---------|-----------|
| `hooks.trainOnSuccess` | **true** | false | Train neural patterns on successful operations |
| `hooks.learnFromFailures` | **true** | false | Learn from failed operations (what NOT to do) |

### Why Learn From Failures?

Negative examples are valuable:
- Prevents repeating mistakes
- Identifies risky patterns
- Improves future decision-making

### Set Hooks Settings

```javascript
mcp__claude-flow__config_set({ key: "hooks.trainOnSuccess", value: true })
mcp__claude-flow__config_set({ key: "hooks.learnFromFailures", value: true })
```

---

## Background Workers

12 background workers handle automated learning tasks.

| Worker | Enabled | Interval | Purpose |
|--------|---------|----------|---------|
| `ultralearn` | **true** | 60s | Deep knowledge acquisition |
| `predict` | **true** | 10min | Predictive preloading |
| `document` | **true** | 60min | Auto-documentation |
| `deepdive` | **true** | 60s | Deep code analysis |
| `refactor` | **true** | 30s | Refactoring suggestions |
| `benchmark` | **true** | 60s | Performance benchmarking |
| `preload` | **true** | 10s | Resource preloading |
| `map` | true | 15min | Codebase mapping |
| `audit` | true | 10min | Security analysis |
| `optimize` | true | 15min | Performance optimization |
| `consolidate` | true | 30min | Memory consolidation (EWC++) |
| `testgaps` | true | 20min | Test coverage analysis |

### Enable Learning Workers

```javascript
// Enable all learning-related workers
mcp__claude-flow__config_set({ key: "daemon.workers.ultralearn.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.predict.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.document.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.deepdive.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.refactor.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.benchmark.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.preload.enabled", value: true })
```

### Daemon Configuration

```javascript
mcp__claude-flow__config_set({ key: "daemon.autoStart", value: true })
mcp__claude-flow__config_set({ key: "daemon.maxConcurrent", value: 16 })
mcp__claude-flow__config_set({ key: "daemon.resourceThresholds.maxCpuLoad", value: 24 })
```

---

## Complete Configuration Reference

### All 62 Settings

```json
{
  "daemon.autoStart": true,
  "daemon.maxConcurrent": 16,
  "daemon.resourceThresholds.maxCpuLoad": 24,
  "daemon.workers.benchmark.enabled": true,
  "daemon.workers.deepdive.enabled": true,
  "daemon.workers.document.enabled": true,
  "daemon.workers.predict.enabled": true,
  "daemon.workers.preload.enabled": true,
  "daemon.workers.refactor.enabled": true,
  "daemon.workers.ultralearn.enabled": true,
  "embeddings.cacheSize": 8192,
  "ewc.enabled": true,
  "ewc.fisherSamples": 100,
  "ewc.gamma": 0.95,
  "ewc.lambda": 0.4,
  "hnsw.efConstruction": 400,
  "hnsw.efSearch": 200,
  "hnsw.m": 48,
  "hnsw.maxElements": 100000000,
  "hooks.learnFromFailures": true,
  "hooks.trainOnSuccess": true,
  "learning.autoLearn": true,
  "learning.batchSize": 32,
  "learning.learningRate": 0.01,
  "learning.patternThreshold": 0.85,
  "logging.format": "json",
  "logging.level": "info",
  "lora.alpha": 32,
  "lora.dropout": 0.05,
  "lora.enabled": true,
  "lora.rank": 16,
  "memory.maxEntries": 50000000,
  "memory.persistInterval": 30000,
  "moe.adaptiveRouting": true,
  "moe.enabled": true,
  "moe.experts": 8,
  "moe.topK": 2,
  "neural.cognitiveModels": 27,
  "neural.enabled": true,
  "neural.epochs": 100,
  "neural.flashAttention": true,
  "patterns.autoStore": true,
  "patterns.clusteringEnabled": true,
  "patterns.enabled": true,
  "performance.memoryLimit": 68719476736,
  "performance.workerThreads": 16,
  "reasoningBank.enabled": true,
  "reasoningBank.verdictTracking": true,
  "security.pathValidation": true,
  "security.sandboxEnabled": true,
  "session.autoSave": true,
  "session.saveInterval": 300000,
  "sona.enabled": true,
  "swarm.autoScale": true,
  "swarm.autoScale.maxAgents": 100,
  "swarm.autoScale.minAgents": 4,
  "swarm.maxAgents": 100,
  "swarm.topology": "hierarchical-mesh",
  "trajectory.autoRecord": true,
  "trajectory.bufferSize": 10000,
  "trajectory.enabled": true,
  "trajectory.maxSteps": 50
}
```

### Quick Copy-Paste: Enable All Learning

```javascript
// SONA
mcp__claude-flow__config_set({ key: "sona.enabled", value: true })

// EWC++
mcp__claude-flow__config_set({ key: "ewc.enabled", value: true })
mcp__claude-flow__config_set({ key: "ewc.lambda", value: 0.4 })
mcp__claude-flow__config_set({ key: "ewc.gamma", value: 0.95 })
mcp__claude-flow__config_set({ key: "ewc.fisherSamples", value: 100 })

// LoRA
mcp__claude-flow__config_set({ key: "lora.enabled", value: true })
mcp__claude-flow__config_set({ key: "lora.rank", value: 16 })
mcp__claude-flow__config_set({ key: "lora.alpha", value: 32 })
mcp__claude-flow__config_set({ key: "lora.dropout", value: 0.05 })

// MoE
mcp__claude-flow__config_set({ key: "moe.enabled", value: true })
mcp__claude-flow__config_set({ key: "moe.adaptiveRouting", value: true })
mcp__claude-flow__config_set({ key: "moe.experts", value: 8 })
mcp__claude-flow__config_set({ key: "moe.topK", value: 2 })

// Trajectories
mcp__claude-flow__config_set({ key: "trajectory.enabled", value: true })
mcp__claude-flow__config_set({ key: "trajectory.autoRecord", value: true })
mcp__claude-flow__config_set({ key: "trajectory.bufferSize", value: 10000 })
mcp__claude-flow__config_set({ key: "trajectory.maxSteps", value: 50 })

// Patterns
mcp__claude-flow__config_set({ key: "patterns.enabled", value: true })
mcp__claude-flow__config_set({ key: "patterns.autoStore", value: true })
mcp__claude-flow__config_set({ key: "patterns.clusteringEnabled", value: true })

// ReasoningBank
mcp__claude-flow__config_set({ key: "reasoningBank.enabled", value: true })
mcp__claude-flow__config_set({ key: "reasoningBank.verdictTracking", value: true })

// Learning
mcp__claude-flow__config_set({ key: "learning.autoLearn", value: true })
mcp__claude-flow__config_set({ key: "learning.learningRate", value: 0.01 })
mcp__claude-flow__config_set({ key: "learning.batchSize", value: 32 })
mcp__claude-flow__config_set({ key: "learning.patternThreshold", value: 0.85 })

// Hooks
mcp__claude-flow__config_set({ key: "hooks.trainOnSuccess", value: true })
mcp__claude-flow__config_set({ key: "hooks.learnFromFailures", value: true })

// Workers
mcp__claude-flow__config_set({ key: "daemon.autoStart", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.ultralearn.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.predict.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.document.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.deepdive.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.refactor.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.benchmark.enabled", value: true })
mcp__claude-flow__config_set({ key: "daemon.workers.preload.enabled", value: true })

// Restart to apply
mcp__claude-flow__system_reset({ component: "all", confirm: true })
```

---

## Troubleshooting

### Learning Not Working?

1. **Check intelligence status:**
   ```javascript
   mcp__claude-flow__hooks_intelligence({ showStatus: true })
   ```

2. **Verify all components are "active":**
   - sona: active
   - moe: active
   - ewc: active
   - lora: active
   - hnsw: active

3. **Check for errors in logs:**
   ```bash
   npx @claude-flow/cli@latest logs --tail 100
   ```

### Patterns Not Being Stored?

- Ensure `patterns.autoStore: true`
- Check `learning.patternThreshold` - patterns below this confidence are ignored
- Verify `hooks.trainOnSuccess: true`

### Memory Issues?

- Reduce `trajectory.bufferSize` if RAM constrained
- Lower `hnsw.m` and `hnsw.maxElements`
- Decrease `embeddings.cacheSize`

---

## References

- [Claude-Flow GitHub](https://github.com/ruvnet/claude-flow)
- [Neural Networks Wiki](https://github.com/ruvnet/claude-flow/wiki/Neural-Networks)
- [Hive-Mind Intelligence Wiki](https://github.com/ruvnet/claude-flow/wiki/Hive-Mind-Intelligence)
- [Self-Improving AI Issue #419](https://github.com/ruvnet/claude-flow/issues/419)
- [EWC Paper](https://arxiv.org/abs/1612.00796)
- [LoRA Paper](https://arxiv.org/abs/2106.09685)

---

*Configuration applied: 2026-01-21*
*Hardware: AMD Ryzen 9 7950X3D, 187GB RAM*
*Total settings: 62*
