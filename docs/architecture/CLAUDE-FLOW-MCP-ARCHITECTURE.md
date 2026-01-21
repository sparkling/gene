# Claude-Flow MCP Server Architecture

## Overview

Claude-Flow is a multi-agent orchestration system that exposes **199 MCP tools** across **23 categories** via a JSON-RPC 2.0 server. It supports multiple transport layers (stdio, HTTP, WebSocket) and integrates with Claude Code through the MCP protocol.

## Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                     Claude Code (Client)                         │
│                                                                  │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐      │
│  │ Task Tool    │    │ MCP Tools    │    │ Bash Tool    │      │
│  │ (Agents)     │    │ (Direct)     │    │ (CLI)        │      │
│  └──────┬───────┘    └──────┬───────┘    └──────┬───────┘      │
└─────────┼──────────────────┼──────────────────┼────────────────┘
          │                   │                   │
          │                   ▼                   │
          │         ┌─────────────────┐           │
          │         │  MCP Protocol   │           │
          │         │  (JSON-RPC 2.0) │           │
          │         │  via stdio      │           │
          │         └────────┬────────┘           │
          │                  │                    │
          ▼                  ▼                    ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Claude-Flow Server                            │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                    Entry Point (bin/cli.js)               │   │
│  │                                                           │   │
│  │   MCP Mode Detection:                                     │   │
│  │   const isMCPMode = !process.stdin.isTTY &&              │   │
│  │                     process.argv.length === 2;            │   │
│  └─────────────────────────┬────────────────────────────────┘   │
│                            │                                     │
│            ┌───────────────┴───────────────┐                    │
│            ▼                               ▼                    │
│  ┌─────────────────┐             ┌─────────────────┐           │
│  │   MCP Server    │             │   CLI Mode      │           │
│  │   (mcp-server)  │             │   (Interactive) │           │
│  └────────┬────────┘             └────────┬────────┘           │
│           │                               │                     │
│           └───────────────┬───────────────┘                    │
│                           │                                     │
│                           ▼                                     │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                    Tool Registry                          │   │
│  │               (TOOL_REGISTRY Map)                         │   │
│  │                                                           │   │
│  │   - O(1) lookup via Map                                   │   │
│  │   - Category/tag indices                                  │   │
│  │   - <10ms registration                                    │   │
│  │   - Execution tracking                                    │   │
│  └─────────────────────────┬────────────────────────────────┘   │
│                            │                                     │
│                            ▼                                     │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                  MCP Tool Modules (26 files)              │   │
│  │                                                           │   │
│  │  agent-tools    swarm-tools    memory-tools   config-tools│   │
│  │  hooks-tools    task-tools     session-tools  hive-mind   │   │
│  │  workflow-tools analyze-tools  progress-tools embeddings  │   │
│  │  claims-tools   security-tools transfer-tools system-tools│   │
│  │  terminal-tools neural-tools   performance    github-tools│   │
│  │  daa-tools      coordination   browser-tools              │   │
│  └──────────────────────────────────────────────────────────┘   │
│                                                                  │
│  ┌──────────────────────────────────────────────────────────┐   │
│  │                    Storage Layer                          │   │
│  │                                                           │   │
│  │  .claude-flow/memory/store.json    (MCP memory)          │   │
│  │  .claude-flow/agents.json          (Agent state)         │   │
│  │  .claude-flow/hive-mind/state.json (Swarm state)         │   │
│  │  .claude-flow/daemon-state.json    (Worker state)        │   │
│  │  .swarm/memory.db                  (CLI SQLite)          │   │
│  └──────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
```

## MCP Mode Detection

The CLI automatically enters MCP server mode when:

```javascript
// bin/cli.js
const isMCPMode = !process.stdin.isTTY && process.argv.length === 2;
```

This enables seamless integration with Claude Code:
```bash
# Claude Code calls this - stdin is piped, no args
claude mcp add claude-flow -- npx -y @claude-flow/cli@latest
```

## Tool Naming Convention

Tools are named using a `category/action` pattern:

| Internal Name | Claude Code Sees |
|---------------|------------------|
| `agent/spawn` | `mcp__claude-flow__agent_spawn` |
| `memory/search` | `mcp__claude-flow__memory_search` |
| `hooks/pre-task` | `mcp__claude-flow__hooks_pre-task` |

## Tool Categories (23 total, 199 tools)

| Category | Tools | Key Functions |
|----------|-------|---------------|
| **hooks** | 38 | Pre/post task, routing, intelligence, workers |
| **browser** | 23 | Page navigation, DOM interaction, screenshots |
| **claims** | 12 | Issue claiming, handoff, load balancing |
| **transfer** | 11 | Pattern store, plugins, PII detection |
| **hive-mind** | 9 | Consensus, broadcast, shared memory |
| **workflow** | 9 | Create, execute, templates |
| **daa** | 8 | Decentralized autonomous agents |
| **agent** | 7 | Spawn, terminate, status, pool |
| **coordination** | 7 | Topology, load balance, sync |
| **embeddings** | 7 | Generate, search, hyperbolic |
| **config** | 6 | Get, set, import, export |
| **memory** | 6 | Store, retrieve, search, stats |
| **neural** | 6 | Train, predict, patterns |
| **performance** | 6 | Benchmark, profile, metrics |
| **task** | 6 | Create, update, complete |
| **analyze** | 6 | Diff risk, reviewers |
| **aidefence** | 6 | Scan, analyze, PII |
| **session** | 5 | Save, restore, list |
| **system** | 5 | Status, health, metrics |
| **terminal** | 5 | Create, execute, history |
| **github** | 5 | Repo analyze, PR, issues |
| **swarm** | 4 | Init, status, shutdown |
| **progress** | 4 | V3 implementation tracking |

## JSON-RPC 2.0 Protocol

### Methods

| Method | Description |
|--------|-------------|
| `initialize` | Returns server info, capabilities, protocol version |
| `tools/list` | Returns all available tools with JSON schemas |
| `tools/call` | Executes a tool with parameters |
| `ping` | Health check |

### Example Message Flow

```json
// Client -> Server (tools/call)
{
  "jsonrpc": "2.0",
  "id": 1,
  "method": "tools/call",
  "params": {
    "name": "memory/store",
    "arguments": {
      "key": "my-key",
      "value": "my-value"
    }
  }
}

// Server -> Client (response)
{
  "jsonrpc": "2.0",
  "id": 1,
  "result": {
    "content": [{
      "type": "text",
      "text": "{\"success\": true, \"key\": \"my-key\"}"
    }]
  }
}
```

## Transport Layers

| Transport | Use Case | Configuration |
|-----------|----------|---------------|
| **stdio** | Claude Code (default) | Automatic via MCP mode detection |
| **HTTP** | REST API clients | `--transport http --port 3000` |
| **WebSocket** | Bidirectional real-time | `--transport websocket` |
| **in-process** | Direct function calls | Programmatic usage |

## Storage Architecture

Claude-Flow uses **multiple storage backends**:

### 1. MCP Memory Store (JSON file)
- **Location**: `.claude-flow/memory/store.json`
- **Purpose**: Simple key-value for MCP tools
- **Schema**: `{ "entries": {...}, "version": "3.0.0" }`

### 2. CLI Memory Store (SQLite)
- **Location**: `.swarm/memory.db`
- **Purpose**: Vector search, HNSW indexing
- **Features**: 384-dim embeddings, semantic search

### 3. State Files
- `.claude-flow/agents.json` - Agent registry
- `.claude-flow/hive-mind/state.json` - Swarm coordination
- `.claude-flow/daemon-state.json` - Background workers

## Model Routing (ADR-026)

3-tier intelligent routing for cost optimization:

| Tier | Handler | Latency | Cost | Use Case |
|------|---------|---------|------|----------|
| 1 | Agent Booster | <1ms | $0 | Code transforms (var→const, add-types) |
| 2 | Haiku | ~500ms | $0.0002 | Simple tasks, bug fixes |
| 3 | Sonnet/Opus | 2-5s | $0.003-$0.015 | Architecture, security |

```javascript
// Model capabilities
const MODEL_CAPABILITIES = {
    haiku: { maxComplexity: 0.4, costMultiplier: 0.04 },
    sonnet: { maxComplexity: 0.7, costMultiplier: 0.2 },
    opus: { maxComplexity: 1.0, costMultiplier: 1.0 },
};
```

## Background Workers (12)

| Worker | Interval | Priority | Purpose |
|--------|----------|----------|---------|
| audit | 10min | critical | Security analysis |
| optimize | 15min | high | Performance optimization |
| map | 15min | normal | Codebase mapping |
| testgaps | 20min | normal | Test coverage |
| consolidate | 30min | low | Memory consolidation |
| ultralearn | - | normal | Deep learning |
| deepdive | - | normal | Code analysis |
| document | - | normal | Auto-documentation |
| refactor | - | normal | Suggestions |
| benchmark | - | normal | Performance |
| predict | - | normal | Preloading |
| preload | - | low | Resource preloading |

## Important Limitations

### GitHub Tools
```javascript
/**
 * ⚠️ IMPORTANT: These tools provide LOCAL STATE MANAGEMENT only.
 * - NO actual GitHub API calls are made
 * - Data is stored locally for workflow coordination
 * - For real GitHub operations, use `gh` CLI or GitHub MCP server
 */
```

### MCP vs CLI Storage
- **MCP tools** use `.claude-flow/memory/store.json` (JSON)
- **CLI tools** use `.swarm/memory.db` (SQLite)
- These are **separate** storage systems

## Configuration Files

| File | Purpose |
|------|---------|
| `.mcp.json` | MCP server registration for Claude Code |
| `.claude/settings.json` | Claude Code project settings |
| `.claude-flow/config.yaml` | Main claude-flow configuration |

## See Also

- [ADR-005: Three-Database Architecture](../adr/ADR-005-THREE-DATABASE-ARCHITECTURE.md)
- [ADR-026: Intelligent Model Routing](../adr/ADR-026-INTELLIGENT-MODEL-ROUTING.md)
- [Tool Inventory](./CLAUDE-FLOW-TOOL-INVENTORY.md)
