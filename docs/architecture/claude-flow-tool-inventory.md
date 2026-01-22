# Claude-Flow MCP Tool Inventory

**Total: 199 tools across 23 categories**

## Quick Reference

| Category | Count | Primary Use |
|----------|-------|-------------|
| hooks | 38 | Task routing, learning, workers |
| browser | 23 | Web automation |
| claims | 12 | Issue management |
| transfer | 11 | Pattern/plugin store |
| hive-mind | 9 | Swarm consensus |
| workflow | 9 | Workflow automation |
| daa | 8 | Autonomous agents |
| agent | 7 | Agent lifecycle |
| coordination | 7 | Swarm topology |
| embeddings | 7 | Vector operations |
| config | 6 | Configuration |
| memory | 6 | Key-value storage |
| neural | 6 | ML operations |
| performance | 6 | Benchmarking |
| task | 6 | Task management |
| analyze | 6 | Git diff analysis |
| aidefence | 6 | Security scanning |
| session | 5 | Session state |
| system | 5 | System status |
| terminal | 5 | Terminal sessions |
| github | 5 | GitHub state (local) |
| swarm | 4 | Swarm lifecycle |
| progress | 4 | V3 tracking |

---

## Complete Tool Listing

### Agent Tools (7)

| Tool | Description |
|------|-------------|
| `agent/spawn` | Spawn agent with intelligent model selection |
| `agent/terminate` | Terminate an agent |
| `agent/status` | Get agent status |
| `agent/list` | List all agents |
| `agent/pool` | Manage agent pool (scale, drain, fill) |
| `agent/health` | Check agent health |
| `agent/update` | Update agent config |

### Swarm Tools (4)

| Tool | Description |
|------|-------------|
| `swarm/init` | Initialize swarm with topology |
| `swarm/status` | Get swarm status |
| `swarm/shutdown` | Shutdown swarm |
| `swarm/health` | Check swarm health |

### Memory Tools (6)

| Tool | Description |
|------|-------------|
| `memory/store` | Store value (persisted to disk) |
| `memory/retrieve` | Retrieve value by key |
| `memory/search` | Semantic search (HNSW-indexed) |
| `memory/delete` | Delete entry |
| `memory/list` | List all entries |
| `memory/stats` | Storage statistics |

### Hooks Tools (38)

| Tool | Description |
|------|-------------|
| `hooks/pre-edit` | Context before file edit |
| `hooks/post-edit` | Record edit outcome |
| `hooks/pre-command` | Risk assessment before command |
| `hooks/post-command` | Record command outcome |
| `hooks/pre-task` | Task start, agent suggestions |
| `hooks/post-task` | Task completion for learning |
| `hooks/route` | Route task to optimal agent |
| `hooks/explain` | Explain routing decision |
| `hooks/metrics` | Learning metrics dashboard |
| `hooks/list` | List registered hooks |
| `hooks/pretrain` | Bootstrap intelligence from repo |
| `hooks/build-agents` | Generate agent configs |
| `hooks/transfer` | Transfer patterns between projects |
| `hooks/session-start` | Initialize session |
| `hooks/session-end` | End session, persist state |
| `hooks/session-restore` | Restore previous session |
| `hooks/notify` | Cross-agent notification |
| `hooks/init` | Initialize hooks in project |
| `hooks/intelligence` | RuVector system status |
| `hooks/intelligence-reset` | Reset learning state |
| `hooks/intelligence/trajectory-start` | Begin SONA trajectory |
| `hooks/intelligence/trajectory-step` | Record trajectory step |
| `hooks/intelligence/trajectory-end` | End trajectory, trigger learning |
| `hooks/intelligence/pattern-store` | Store pattern (HNSW) |
| `hooks/intelligence/pattern-search` | Search patterns |
| `hooks/intelligence/stats` | Intelligence statistics |
| `hooks/intelligence/learn` | Force learning cycle |
| `hooks/intelligence/attention` | Attention-weighted similarity |
| `hooks/worker-list` | List 12 background workers |
| `hooks/worker-dispatch` | Dispatch worker |
| `hooks/worker-status` | Worker status |
| `hooks/worker-detect` | Detect triggers from prompt |
| `hooks/worker-cancel` | Cancel worker |
| `hooks/model-route` | Route to haiku/sonnet/opus |
| `hooks/model-outcome` | Record model outcome |
| `hooks/model-stats` | Model routing statistics |

### Hive-Mind Tools (9)

| Tool | Description |
|------|-------------|
| `hive-mind/spawn` | Spawn workers and join hive |
| `hive-mind/init` | Initialize collective |
| `hive-mind/status` | Hive status |
| `hive-mind/join` | Join agent to hive |
| `hive-mind/leave` | Remove agent |
| `hive-mind/consensus` | Propose/vote on consensus |
| `hive-mind/broadcast` | Broadcast to all workers |
| `hive-mind/shutdown` | Shutdown hive |
| `hive-mind/memory` | Shared memory access |

### Browser Tools (23)

| Tool | Description |
|------|-------------|
| `browser/open` | Navigate to URL |
| `browser/back` | Navigate back |
| `browser/forward` | Navigate forward |
| `browser/reload` | Reload page |
| `browser/close` | Close session |
| `browser/snapshot` | Accessibility tree |
| `browser/screenshot` | Capture screenshot |
| `browser/click` | Click element |
| `browser/fill` | Fill input |
| `browser/type` | Type with key events |
| `browser/press` | Press key |
| `browser/hover` | Hover element |
| `browser/select` | Select dropdown |
| `browser/check` | Check checkbox |
| `browser/uncheck` | Uncheck checkbox |
| `browser/scroll` | Scroll page |
| `browser/get-text` | Get element text |
| `browser/get-value` | Get input value |
| `browser/get-title` | Get page title |
| `browser/get-url` | Get current URL |
| `browser/wait` | Wait for condition |
| `browser/eval` | Execute JavaScript |
| `browser/session-list` | List sessions |

### Claims Tools (12)

| Tool | Description |
|------|-------------|
| `claims/claim` | Claim issue for work |
| `claims/release` | Release claim |
| `claims/handoff` | Request handoff |
| `claims/accept-handoff` | Accept handoff |
| `claims/status` | Update claim status |
| `claims/list` | List all claims |
| `claims/mark-stealable` | Mark as stealable |
| `claims/steal` | Steal issue |
| `claims/stealable` | List stealable issues |
| `claims/load` | Agent load info |
| `claims/board` | Visual board view |
| `claims/rebalance` | Load rebalancing |

### Config Tools (6)

| Tool | Description |
|------|-------------|
| `config/get` | Get config value |
| `config/set` | Set config value |
| `config/list` | List config |
| `config/reset` | Reset to defaults |
| `config/export` | Export to JSON |
| `config/import` | Import from JSON |

### Task Tools (6)

| Tool | Description |
|------|-------------|
| `task/create` | Create task |
| `task/status` | Task status |
| `task/list` | List tasks |
| `task/complete` | Mark complete |
| `task/update` | Update task |
| `task/cancel` | Cancel task |

### Session Tools (5)

| Tool | Description |
|------|-------------|
| `session/save` | Save session state |
| `session/restore` | Restore session |
| `session/list` | List sessions |
| `session/delete` | Delete session |
| `session/info` | Session info |

### Workflow Tools (9)

| Tool | Description |
|------|-------------|
| `workflow/create` | Create workflow |
| `workflow/execute` | Execute workflow |
| `workflow/status` | Workflow status |
| `workflow/list` | List workflows |
| `workflow/pause` | Pause workflow |
| `workflow/resume` | Resume workflow |
| `workflow/cancel` | Cancel workflow |
| `workflow/delete` | Delete workflow |
| `workflow/template` | Save/create templates |

### Analyze Tools (6)

| Tool | Description |
|------|-------------|
| `analyze/diff` | Git diff risk assessment |
| `analyze/diff-risk` | Quick risk assessment |
| `analyze/diff-classify` | Classify change type |
| `analyze/diff-reviewers` | Suggest reviewers |
| `analyze/file-risk` | File change risk |
| `analyze/diff-stats` | Diff statistics |

### Embeddings Tools (7)

| Tool | Description |
|------|-------------|
| `embeddings/init` | Initialize ONNX subsystem |
| `embeddings/generate` | Generate embeddings |
| `embeddings/compare` | Compare similarity |
| `embeddings/search` | Semantic search |
| `embeddings/neural` | Neural substrate ops |
| `embeddings/hyperbolic` | Poincare ball ops |
| `embeddings/status` | System status |

### Neural Tools (6)

| Tool | Description |
|------|-------------|
| `neural/train` | Train model |
| `neural/predict` | Make predictions |
| `neural/patterns` | Manage patterns |
| `neural/compress` | Compress model |
| `neural/status` | Neural status |
| `neural/optimize` | Optimize model |

### Performance Tools (6)

| Tool | Description |
|------|-------------|
| `performance/report` | Generate report |
| `performance/bottleneck` | Detect bottlenecks |
| `performance/benchmark` | Run benchmarks |
| `performance/profile` | Profile component |
| `performance/optimize` | Apply optimizations |
| `performance/metrics` | Get metrics |

### System Tools (5)

| Tool | Description |
|------|-------------|
| `system/status` | Overall status |
| `system/metrics` | System metrics |
| `system/health` | Health check |
| `system/info` | System info |
| `system/reset` | Reset state |

### Terminal Tools (5)

| Tool | Description |
|------|-------------|
| `terminal/create` | Create session |
| `terminal/execute` | Execute command |
| `terminal/list` | List sessions |
| `terminal/close` | Close session |
| `terminal/history` | Command history |

### GitHub Tools (5)

> **Note**: LOCAL STATE ONLY - no actual GitHub API calls

| Tool | Description |
|------|-------------|
| `github/repo_analyze` | Analyze repository |
| `github/pr_manage` | Manage PRs |
| `github/issue_track` | Track issues |
| `github/workflow` | Manage workflows |
| `github/metrics` | Repository metrics |

### AIDefence Tools (6)

| Tool | Description |
|------|-------------|
| `aidefence_scan` | Scan for threats |
| `aidefence_analyze` | Deep threat analysis |
| `aidefence_stats` | Detection stats |
| `aidefence_learn` | Record feedback |
| `aidefence_is_safe` | Quick safety check |
| `aidefence_has_pii` | PII detection |

### Transfer Tools (11)

| Tool | Description |
|------|-------------|
| `transfer/detect-pii` | Detect PII |
| `transfer/ipfs-resolve` | Resolve IPNS name |
| `transfer/store-search` | Search pattern store |
| `transfer/store-info` | Pattern info |
| `transfer/store-download` | Download pattern |
| `transfer/store-featured` | Featured patterns |
| `transfer/store-trending` | Trending patterns |
| `transfer/plugin-search` | Search plugins |
| `transfer/plugin-info` | Plugin info |
| `transfer/plugin-featured` | Featured plugins |
| `transfer/plugin-official` | Official plugins |

### Progress Tools (4)

| Tool | Description |
|------|-------------|
| `progress/check` | V3 progress percentage |
| `progress/sync` | Persist progress metrics |
| `progress/summary` | Human-readable summary |
| `progress/watch` | Watch status |

### Coordination Tools (7)

| Tool | Description |
|------|-------------|
| `coordination/topology` | Configure topology |
| `coordination/load_balance` | Load balancing |
| `coordination/sync` | State sync |
| `coordination/node` | Manage nodes |
| `coordination/consensus` | Consensus protocol |
| `coordination/orchestrate` | Multi-agent coordination |
| `coordination/metrics` | Coordination metrics |

### DAA Tools (8)

| Tool | Description |
|------|-------------|
| `daa/agent_create` | Create autonomous agent |
| `daa/agent_adapt` | Trigger adaptation |
| `daa/workflow_create` | Create workflow |
| `daa/workflow_execute` | Execute workflow |
| `daa/knowledge_share` | Share knowledge |
| `daa/learning_status` | Learning status |
| `daa/cognitive_pattern` | Cognitive patterns |
| `daa/performance_metrics` | Performance metrics |

---

## Tool Registration

Tools are registered via a centralized Map:

```javascript
const TOOL_REGISTRY = new Map();

function registerTools(tools) {
    tools.forEach(tool => TOOL_REGISTRY.set(tool.name, tool));
}
```

## Tool Schema

Each tool follows this interface:

```typescript
interface MCPTool {
  name: string;           // e.g., "agent/spawn"
  description: string;
  category?: string;
  inputSchema: {
    type: 'object';
    properties: Record<string, any>;
    required?: string[];
  };
  handler: (input: any) => Promise<any>;
}
```

## Configuration

### Enable Specific Tools

```bash
npx @claude-flow/cli@latest mcp start --tools agent/spawn,memory/search
```

### Runtime Toggle

```bash
npx @claude-flow/cli@latest mcp toggle --enable agent/spawn
npx @claude-flow/cli@latest mcp toggle --disable browser/eval
```
