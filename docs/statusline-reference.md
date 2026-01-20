# Claude Code Statusline Reference

Complete reference for statusline configuration, data sources, and design decisions.

## Table of Contents

1. [Claude Code Native Statusline](#claude-code-native-statusline)
2. [Configuration](#configuration)
3. [Data Available via stdin](#data-available-via-stdin)
4. [Claude-Flow Data Sources](#claude-flow-data-sources)
5. [Third-Party Tools](#third-party-tools)
6. [All Available Data Points](#all-available-data-points)
7. [Mode Design & Detection](#mode-design--detection)
8. [Implementation Notes](#implementation-notes)

---

## Claude Code Native Statusline

### Overview
- Displayed at bottom of Claude Code terminal
- Updates when conversation messages change (max every 300ms)
- Configurable via `.claude/settings.json` or `/statusline` command
- Script receives JSON via stdin, outputs first line of stdout

### `/statusline` Command
```
/statusline [optional instructions]
```
Interactive setup wizard for configuring status line script.

---

## Configuration

### settings.json
```json
{
  "statusLine": {
    "type": "command",
    "command": "bash .claude/statusline.sh",
    "padding": 0,
    "enabled": true,
    "refreshMs": 5000
  }
}
```

| Option | Description |
|--------|-------------|
| `type` | Set to `"command"` for custom script |
| `command` | Path to executable script |
| `padding` | Set to 0 for edge-to-edge display |
| `enabled` | Toggle statusline on/off |
| `refreshMs` | Refresh interval in milliseconds |

---

## Data Available via stdin

Claude Code sends this JSON to the statusline script:

```json
{
  "hook_event_name": "Status",
  "session_id": "abc123...",
  "transcript_path": "/path/to/transcript.json",
  "cwd": "/current/working/directory",
  "model": {
    "id": "claude-opus-4-1",
    "display_name": "Opus"
  },
  "workspace": {
    "current_dir": "/current/working/directory",
    "project_dir": "/original/project/directory"
  },
  "version": "1.0.80",
  "output_style": {
    "name": "default"
  },
  "cost": {
    "total_cost_usd": 0.01234,
    "total_duration_ms": 45000,
    "total_api_duration_ms": 2300,
    "total_lines_added": 156,
    "total_lines_removed": 23
  },
  "context_window": {
    "total_input_tokens": 15234,
    "total_output_tokens": 4521,
    "context_window_size": 200000,
    "used_percentage": 42.5,
    "remaining_percentage": 57.5,
    "current_usage": {
      "input_tokens": 8500,
      "output_tokens": 1200,
      "cache_creation_input_tokens": 5000,
      "cache_read_input_tokens": 2000
    }
  }
}
```

### Field Reference

| Category | Field | Description |
|----------|-------|-------------|
| **Model** | `model.display_name` | Human-readable name (e.g., "Opus 4.5") |
| | `model.id` | Full model identifier |
| **Workspace** | `workspace.current_dir` | Current working directory |
| | `workspace.project_dir` | Original project directory |
| **Context** | `context_window.used_percentage` | Context usage 0-100 |
| | `context_window.remaining_percentage` | Context remaining 0-100 |
| | `context_window.context_window_size` | Total context size |
| | `context_window.total_input_tokens` | Cumulative input tokens |
| | `context_window.total_output_tokens` | Cumulative output tokens |
| **Cost** | `cost.total_cost_usd` | Session cost in USD |
| | `cost.total_duration_ms` | Total session duration |
| | `cost.total_api_duration_ms` | Time in API calls |
| **Code** | `cost.total_lines_added` | Lines added this session |
| | `cost.total_lines_removed` | Lines removed this session |
| **Session** | `session_id` | Unique session identifier |
| | `transcript_path` | Path to conversation transcript |
| | `version` | Claude Code version |

---

## Claude-Flow Data Sources

### CLI Statusline Hook
```bash
npx @claude-flow/cli@latest hooks statusline --json
```

Returns:
```json
{
  "user": {
    "name": "Claude",
    "gitBranch": "main",
    "modelName": "Opus 4.5"
  },
  "v3Progress": {
    "domainsCompleted": 3,
    "totalDomains": 5,
    "dddProgress": 60,
    "patternsLearned": 154,
    "sessionsCompleted": 15
  },
  "security": {
    "status": "PENDING",
    "cvesFixed": 0,
    "totalCves": 3
  },
  "swarm": {
    "activeAgents": 1,
    "maxAgents": 15,
    "coordinationActive": true
  },
  "system": {
    "memoryMB": 18,
    "contextPct": 75,
    "intelligencePct": 15,
    "subAgents": 0
  },
  "timestamp": "2026-01-20T09:55:17.487Z"
}
```

### File: `.ruvector/intelligence.json`
```json
{
  "stats": {
    "total_patterns": 2,
    "total_memories": 477,
    "total_trajectories": 293,
    "total_errors": 0,
    "session_count": 7
  },
  "agents": [],
  "edges": [],
  "file_sequences": [],
  "memories": [],
  "patterns": [],
  "trajectories": []
}
```

### File: `.claude-flow/daemon-state.json`
```json
{
  "running": true,
  "startedAt": "2026-01-20T09:33:44.611Z",
  "workers": {
    "map": {
      "runCount": 210,
      "successCount": 210,
      "failureCount": 0,
      "averageDurationMs": 0.37,
      "lastRun": "timestamp",
      "nextRun": "timestamp",
      "isRunning": false
    }
    // ... other workers: audit, optimize, consolidate, testgaps, predict, etc.
  }
}
```

### File: `.claude-flow/metrics/performance.json`
```json
{
  "flashAttention": {
    "enabled": true,
    "speedup": "1.0x"
  },
  "hnsw": {
    "searchTimeMs": 3,
    "indexSize": 12
  },
  "routing": {
    "costSavings": "51.6%"
  }
}
```

### File: `.claude-flow/security/audit-status.json`
```json
{
  "status": "CLEAN",
  "cvesFixed": 0,
  "totalCves": 0,
  "lastScan": "2026-01-20T02:00:00Z",
  "issues": {
    "critical": 0,
    "high": 0,
    "medium": 0,
    "low": 0
  }
}
```

### Database Files
- `data/user.db` - User graph data (768 dims, WAL mode)
- `data/operational.db` - Claude-flow operational data (384 dims)

Estimate vectors: `file_size_bytes / (dimensions * 2)`

---

## Third-Party Tools

### ccstatusline (by @sirmalloc)
- **GitHub**: https://github.com/sirmalloc/ccstatusline
- **Install**: `bunx ccstatusline@latest`
- **Features**:
  - Up to 3 status lines
  - Powerline mode with themes
  - Interactive TUI configuration
  - Widgets: model, git, tokens, context %, cost, elapsed time, CWD

### claude-code-statusline (by @rz1989s)
- **GitHub**: https://github.com/rz1989s/claude-code-statusline
- **Features**:
  - 18 atomic components across 1-9 lines
  - MCP server health monitoring
  - Cost tracking (daily/weekly/monthly)
  - Burn rate, cache efficiency
  - 18+ themes including Catppuccin
  - Config.toml with 227 settings

### Comparison
| Feature | ccstatusline | claude-code-statusline | Our Custom |
|---------|--------------|------------------------|------------|
| Lines | 3 | 1-9 | Adaptive |
| Claude-flow aware | No | No | Yes |
| Swarm status | No | No | Yes |
| Vector DB stats | No | No | Yes |
| Learning metrics | No | No | Yes |
| MCP health | No | Yes | No |
| Cost tracking | Yes | Yes | Via stdin |

---

## All Available Data Points

**Total: 170+ unique data points across 13 sources**

### Summary

| Source | File/Command | Count |
|--------|--------------|-------|
| Claude Code stdin | (native JSON) | 24 |
| RuVector Intelligence | `.ruvector/intelligence.json` | 15+ |
| Daemon State | `.claude-flow/daemon-state.json` | 55+ |
| V3 Progress | `.claude-flow/metrics/v3-progress.json` | 9 |
| Performance Metrics | `.claude-flow/metrics/performance.json` | 5 |
| Consolidation Metrics | `.claude-flow/metrics/consolidation.json` | 4 |
| Security Audit | `.claude-flow/security/audit-status.json` | 8 |
| Embeddings Config | `.claude-flow/embeddings.json` | 12 |
| Claude-Flow Config | `.claude-flow/config.json` | 13 |
| Codebase Map | `.claude-flow/metrics/codebase-map.json` | 7 |
| Git | (commands) | 6 |
| Database Files | `data/*.db` | 6 |
| System | (commands) | 6 |

---

### By Purpose

Data points organized by display use case:

#### Identity (7 points)
What/where am I working?

| Data Point | Source | Example |
|------------|--------|---------|
| Model display name | `stdin → model.display_name` | `"Opus 4.5"` |
| Model ID | `stdin → model.id` | `"claude-opus-4-5-20251101"` |
| Current directory | `stdin → workspace.current_dir` | `"/home/claude/src/gene"` |
| Project directory | `stdin → workspace.project_dir` | `"/home/claude/src/gene"` |
| Git branch | `git branch --show-current` | `"main"` |
| Session ID | `stdin → session_id` | `"abc123..."` |
| Hostname | `hostname` | `"dev-server"` |

#### Cost & Usage (12 points)
How much am I spending/using?

| Data Point | Source | Example |
|------------|--------|---------|
| Session cost USD | `stdin → cost.total_cost_usd` | `$0.0123` |
| Context used % | `stdin → context_window.used_percentage` | `42.5%` |
| Context remaining % | `stdin → context_window.remaining_percentage` | `57.5%` |
| Context window size | `stdin → context_window.context_window_size` | `200000` |
| Total input tokens | `stdin → context_window.total_input_tokens` | `15234` |
| Total output tokens | `stdin → context_window.total_output_tokens` | `4521` |
| Cache creation tokens | `stdin → ...cache_creation_input_tokens` | `5000` |
| Cache read tokens | `stdin → ...cache_read_input_tokens` | `2000` |
| Session duration | `stdin → cost.total_duration_ms` | `45000ms` |
| API duration | `stdin → cost.total_api_duration_ms` | `2300ms` |
| Lines added | `stdin → cost.total_lines_added` | `156` |
| Lines removed | `stdin → cost.total_lines_removed` | `23` |

#### Learning & Intelligence (20+ points)
What has the system learned?

| Data Point | Source | Example |
|------------|--------|---------|
| Total patterns | `intelligence.json → stats.total_patterns` | `2` |
| Total memories | `intelligence.json → stats.total_memories` | `477` |
| Total trajectories | `intelligence.json → stats.total_trajectories` | `293` |
| Total errors | `intelligence.json → stats.total_errors` | `0` |
| Session count | `intelligence.json → stats.session_count` | `7` |
| Pattern Q-values | `intelligence.json → patterns[].q_value` | `0.799` |
| Pattern visits | `intelligence.json → patterns[].visits` | `160` |
| Patterns learned | `v3-progress.json → patternsLearned` | `140` |
| Sessions completed | `v3-progress.json → sessionsCompleted` | `14` |
| Patterns consolidated | `consolidation.json → patternsConsolidated` | `0` |
| Memory cleaned | `consolidation.json → memoryCleaned` | `0` |
| Duplicates removed | `consolidation.json → duplicatesRemoved` | `0` |
| Neural enabled | `embeddings.json → neural.enabled` | `true` |
| Neural drift threshold | `embeddings.json → neural.driftThreshold` | `0.3` |
| Hyperbolic enabled | `embeddings.json → hyperbolic.enabled` | `true` |

#### Swarm & Agents (15+ points)
What agents are running?

| Data Point | Source | Example |
|------------|--------|---------|
| Active agents | `v3-progress.json → swarm.activeAgents` | `0` |
| Max agents | `v3-progress.json → swarm.maxAgents` | `50` |
| Swarm topology | `config.json → values.swarm.topology` | `"hierarchical-mesh"` |
| Auto scale | `config.json → values.swarm.autoScale` | `true` |
| Daemon running | `daemon-state.json → running` | `true` |
| Daemon started | `daemon-state.json → startedAt` | (timestamp) |
| Worker run counts | `daemon-state.json → workers.*.runCount` | `234` |
| Worker success | `daemon-state.json → workers.*.successCount` | `234` |
| Worker failures | `daemon-state.json → workers.*.failureCount` | `0` |
| Worker is running | `daemon-state.json → workers.*.isRunning` | `false` |
| Max concurrent | `daemon-state.json → config.maxConcurrent` | `2` |

#### Security (8 points)
Are there vulnerabilities?

| Data Point | Source | Example |
|------------|--------|---------|
| Audit status | `audit-status.json → status` | `"CLEAN"` |
| CVEs fixed | `audit-status.json → cvesFixed` | `0` |
| Total CVEs | `audit-status.json → totalCves` | `0` |
| Last scan | `audit-status.json → lastScan` | (timestamp) |
| Critical issues | `audit-status.json → issues.critical` | `0` |
| High issues | `audit-status.json → issues.high` | `0` |
| Medium issues | `audit-status.json → issues.medium` | `0` |
| Low issues | `audit-status.json → issues.low` | `0` |

#### Performance (8 points)
How fast is the system?

| Data Point | Source | Example |
|------------|--------|---------|
| Flash attention enabled | `performance.json → flashAttention.enabled` | `true` |
| Flash attention speedup | `performance.json → flashAttention.speedup` | `"1.0x"` |
| HNSW search time | `performance.json → hnsw.searchTimeMs` | `3ms` |
| HNSW index size | `performance.json → hnsw.indexSize` | `12` |
| Routing cost savings | `performance.json → routing.costSavings` | `"51.6%"` |
| Worker avg duration | `daemon-state.json → workers.*.averageDurationMs` | `0.32ms` |
| Memory persist interval | `config.json → values.memory.persistInterval` | `60000ms` |
| Session save interval | `config.json → values.session.saveInterval` | `300000ms` |

#### Storage & Database (15 points)
What data is persisted?

| Data Point | Source | Example |
|------------|--------|---------|
| User DB size | `du -h data/user.db` | `"1.2M"` |
| User vectors (est) | `bytes / 1536` | `819` |
| Ops DB size | `du -h data/operational.db` | `"512K"` |
| Ops vectors (est) | `bytes / 768` | `682` |
| Memory entries | `v3-progress.json → memory.entries` | `12` |
| Memory size MB | `v3-progress.json → memory.sizeMb` | `18` |
| Memory max entries | `config.json → values.memory.maxEntries` | `1000000` |
| Embedding model | `embeddings.json → model` | `"all-MiniLM-L6-v2"` |
| Embedding dimension | `embeddings.json → dimension` | `384` |
| Cache size | `embeddings.json → cacheSize` | `256` |
| Disk space | `df -h .` | `"45G available"` |

#### Project State (10 points)
What's the project status?

| Data Point | Source | Example |
|------------|--------|---------|
| Domains completed | `v3-progress.json → domains.completed` | `3` |
| Domains total | `v3-progress.json → domains.total` | `5` |
| DDD progress % | `v3-progress.json → ddd.progress` | `60%` |
| Has package.json | `codebase-map.json → structure.hasPackageJson` | `true` |
| Has tsconfig | `codebase-map.json → structure.hasTsConfig` | `false` |
| Has Claude config | `codebase-map.json → structure.hasClaudeConfig` | `true` |
| Has Claude-Flow | `codebase-map.json → structure.hasClaudeFlow` | `true` |
| Uncommitted changes | `git status --porcelain \| wc -l` | `2` |
| Last commit | `git log --oneline -1` | `"e249b3d docs:..."` |
| Stash count | `git stash list \| wc -l` | `0` |

#### System (5 points)
What's the environment?

| Data Point | Source | Example |
|------------|--------|---------|
| Current time | `date` | `"Mon Jan 20..."` |
| System uptime | `uptime` | `"up 5 days"` |
| Memory usage | `free -m` | (stats) |
| Current user | `whoami` | `"claude"` |
| Claude Code version | `stdin → version` | `"1.0.80"` |

---

### By Source

Detailed listing of all data points organized by their source file/command:

### 1. Claude Code stdin (24 points)

Sent automatically to statusline script via stdin JSON.

| Data Point | Source | Example |
|------------|--------|---------|
| Model ID | `model.id` | `"claude-opus-4-5-20251101"` |
| Model display name | `model.display_name` | `"Opus 4.5"` |
| Current directory | `workspace.current_dir` | `"/home/claude/src/gene"` |
| Project directory | `workspace.project_dir` | `"/home/claude/src/gene"` |
| Working directory | `cwd` | `"/home/claude/src/gene"` |
| Session ID | `session_id` | `"abc123..."` |
| Transcript path | `transcript_path` | `"/path/to/transcript.json"` |
| Claude Code version | `version` | `"1.0.80"` |
| Output style | `output_style.name` | `"default"` |
| Context window size | `context_window.context_window_size` | `200000` |
| Context used % | `context_window.used_percentage` | `42.5` |
| Context remaining % | `context_window.remaining_percentage` | `57.5` |
| Total input tokens | `context_window.total_input_tokens` | `15234` |
| Total output tokens | `context_window.total_output_tokens` | `4521` |
| Current input tokens | `context_window.current_usage.input_tokens` | `8500` |
| Current output tokens | `context_window.current_usage.output_tokens` | `1200` |
| Cache creation tokens | `context_window.current_usage.cache_creation_input_tokens` | `5000` |
| Cache read tokens | `context_window.current_usage.cache_read_input_tokens` | `2000` |
| Session cost USD | `cost.total_cost_usd` | `0.01234` |
| Session duration ms | `cost.total_duration_ms` | `45000` |
| API duration ms | `cost.total_api_duration_ms` | `2300` |
| Lines added | `cost.total_lines_added` | `156` |
| Lines removed | `cost.total_lines_removed` | `23` |
| Hook event name | `hook_event_name` | `"Status"` |

---

### 2. RuVector Intelligence (15+ points)

File: `.ruvector/intelligence.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Total patterns | `stats.total_patterns` | `2` |
| Total memories | `stats.total_memories` | `477` |
| Total trajectories | `stats.total_trajectories` | `293` |
| Total errors | `stats.total_errors` | `0` |
| Session count | `stats.session_count` | `7` |
| Last session | `stats.last_session` | `1768871732` |
| Agents count | `agents.length` | (array length) |
| Edges count | `edges.length` | (array length) |
| File sequences count | `file_sequences.length` | (array length) |
| Errors array | `errors.length` | (array length) |

**Pattern entries** (per pattern):
| Field | Source | Example |
|-------|--------|---------|
| State | `patterns[key].state` | `"cmd_shell_general"` |
| Action | `patterns[key].action` | `"success"` |
| Q-value | `patterns[key].q_value` | `0.7999...` |
| Visits | `patterns[key].visits` | `160` |
| Last update | `patterns[key].last_update` | `1768871846` |

**Memory entries** (per memory):
| Field | Source | Example |
|-------|--------|---------|
| ID | `memories[].id` | `"mem_1768860137"` |
| Type | `memories[].memory_type` | `"agent_spawn"` |
| Content | `memories[].content` | (text content) |
| Embedding | `memories[].embedding` | (vector array) |

---

### 3. Daemon State (55+ points)

File: `.claude-flow/daemon-state.json`

**Base fields:**
| Data Point | Source | Example |
|------------|--------|---------|
| Daemon running | `running` | `true` |
| Started at | `startedAt` | `"2026-01-18T16:49:15.823Z"` |
| Saved at | `savedAt` | `"2026-01-20T15:53:08.616Z"` |

**Per worker** (7 workers × 7 fields = 49 points):

Workers: `map`, `audit`, `optimize`, `consolidate`, `testgaps`, `predict`, `document`

| Field | Source | Example |
|-------|--------|---------|
| Run count | `workers.{name}.runCount` | `234` |
| Success count | `workers.{name}.successCount` | `234` |
| Failure count | `workers.{name}.failureCount` | `0` |
| Avg duration ms | `workers.{name}.averageDurationMs` | `0.32` |
| Last run | `workers.{name}.lastRun` | `"2026-01-20T15:53:08.616Z"` |
| Next run | `workers.{name}.nextRun` | `"2026-01-20T16:04:05.680Z"` |
| Is running | `workers.{name}.isRunning` | `false` |

**Config fields:**
| Data Point | Source | Example |
|------------|--------|---------|
| Auto start | `config.autoStart` | `false` |
| Log directory | `config.logDir` | `"/home/.../logs"` |
| State file | `config.stateFile` | `"/home/.../daemon-state.json"` |
| Max concurrent | `config.maxConcurrent` | `2` |
| Worker timeout ms | `config.workerTimeoutMs` | `300000` |
| Max CPU load | `config.resourceThresholds.maxCpuLoad` | `2` |
| Min free memory % | `config.resourceThresholds.minFreeMemoryPercent` | `20` |

---

### 4. V3 Progress (9 points)

File: `.claude-flow/metrics/v3-progress.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Domains completed | `domains.completed` | `3` |
| Domains total | `domains.total` | `5` |
| DDD progress % | `ddd.progress` | `60` |
| Swarm active agents | `swarm.activeAgents` | `0` |
| Swarm max agents | `swarm.maxAgents` | `50` |
| Memory entries | `memory.entries` | `12` |
| Memory size MB | `memory.sizeMb` | `18` |
| Patterns learned | `patternsLearned` | `140` |
| Sessions completed | `sessionsCompleted` | `14` |

---

### 5. Performance Metrics (5 points)

File: `.claude-flow/metrics/performance.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Flash attention enabled | `flashAttention.enabled` | `true` |
| Flash attention speedup | `flashAttention.speedup` | `"1.0x"` |
| HNSW search time ms | `hnsw.searchTimeMs` | `3` |
| HNSW index size | `hnsw.indexSize` | `12` |
| Routing cost savings | `routing.costSavings` | `"51.6%"` |

---

### 6. Consolidation Metrics (4 points)

File: `.claude-flow/metrics/consolidation.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Timestamp | `timestamp` | `"2026-01-20T19:43:34.116Z"` |
| Patterns consolidated | `patternsConsolidated` | `0` |
| Memory cleaned | `memoryCleaned` | `0` |
| Duplicates removed | `duplicatesRemoved` | `0` |

---

### 7. Security Audit (8 points)

File: `.claude-flow/security/audit-status.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Status | `status` | `"CLEAN"` |
| CVEs fixed | `cvesFixed` | `0` |
| Total CVEs | `totalCves` | `0` |
| Last scan | `lastScan` | `"2026-01-20T02:00:00Z"` |
| Critical issues | `issues.critical` | `0` |
| High issues | `issues.high` | `0` |
| Medium issues | `issues.medium` | `0` |
| Low issues | `issues.low` | `0` |

---

### 8. Embeddings Config (12 points)

File: `.claude-flow/embeddings.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Model | `model` | `"all-MiniLM-L6-v2"` |
| Model path | `modelPath` | `"/home/.../models"` |
| Dimension | `dimension` | `384` |
| Cache size | `cacheSize` | `256` |
| Hyperbolic enabled | `hyperbolic.enabled` | `true` |
| Hyperbolic curvature | `hyperbolic.curvature` | `-1` |
| Hyperbolic epsilon | `hyperbolic.epsilon` | `1e-15` |
| Hyperbolic max norm | `hyperbolic.maxNorm` | `0.99999` |
| Neural enabled | `neural.enabled` | `true` |
| Neural drift threshold | `neural.driftThreshold` | `0.3` |
| Neural decay rate | `neural.decayRate` | `0.01` |
| Initialized | `initialized` | `"2026-01-20T01:43:09.406Z"` |

---

### 9. Claude-Flow Config (13 points)

File: `.claude-flow/config.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Swarm topology | `values.swarm.topology` | `"hierarchical-mesh"` |
| Swarm max agents | `values.swarm.maxAgents` | `50` |
| Swarm auto scale | `values.swarm.autoScale` | `true` |
| Memory persist interval | `values.memory.persistInterval` | `60000` |
| Memory max entries | `values.memory.maxEntries` | `1000000` |
| Session auto save | `values.session.autoSave` | `true` |
| Session save interval | `values.session.saveInterval` | `300000` |
| Logging level | `values.logging.level` | `"info"` |
| Logging format | `values.logging.format` | `"json"` |
| Security sandbox enabled | `values.security.sandboxEnabled` | `true` |
| Security path validation | `values.security.pathValidation` | `true` |
| Version | `version` | `"3.0.0"` |
| Updated at | `updatedAt` | `"2026-01-20T01:54:09.593Z"` |

---

### 10. Codebase Map (7 points)

File: `.claude-flow/metrics/codebase-map.json`

| Data Point | Source | Example |
|------------|--------|---------|
| Timestamp | `timestamp` | `"2026-01-20T15:53:08.616Z"` |
| Project root | `projectRoot` | `"/home/claude/src/gene"` |
| Scanned at | `scannedAt` | `1768924388616` |
| Has package.json | `structure.hasPackageJson` | `true` |
| Has tsconfig | `structure.hasTsConfig` | `false` |
| Has Claude config | `structure.hasClaudeConfig` | `true` |
| Has Claude-Flow | `structure.hasClaudeFlow` | `true` |

---

### 11. Git (6 points)

Via shell commands:

| Data Point | Command | Example |
|------------|---------|---------|
| Current branch | `git branch --show-current` | `"main"` |
| Uncommitted changes | `git status --porcelain \| wc -l` | `2` |
| Last commit | `git log --oneline -1` | `"e249b3d docs: expand..."` |
| Stash count | `git stash list \| wc -l` | `0` |
| Short commit hash | `git rev-parse --short HEAD` | `"e249b3d"` |
| Remote URLs | `git remote -v` | (remote list) |

---

### 12. Database Files (6 points)

Files: `data/user.db`, `data/operational.db`

| Data Point | Command/Formula | Example |
|------------|-----------------|---------|
| User DB size | `du -h data/user.db` | `"1.2M"` |
| User DB bytes | `stat -c%s data/user.db` | `1258291` |
| User vectors (est) | `bytes / 1536` | `819` |
| Ops DB size | `du -h data/operational.db` | `"512K"` |
| Ops DB bytes | `stat -c%s data/operational.db` | `524288` |
| Ops vectors (est) | `bytes / 768` | `682` |

---

### 13. System (6 points)

Via shell commands:

| Data Point | Command | Example |
|------------|---------|---------|
| Current time | `date` | `"Mon Jan 20 21:00:00 UTC 2026"` |
| System uptime | `uptime` | `"up 5 days, 3:42"` |
| Memory usage | `free -m` | (memory stats) |
| Disk space | `df -h .` | `"45G available"` |
| Current user | `whoami` | `"claude"` |
| Hostname | `hostname` | `"dev-server"` |

---

## Mode Design & Detection

### Critique of Mode-Based Approach

**Problems identified:**
1. Most "modes" are not mutually exclusive
2. Most "modes" are accumulated data, not transient states
3. Detection can be unreliable (workers finish before refresh)

**Actually transient states:**
- Swarm running (`activeAgents > 1`)
- Alert condition (CVEs > 0 or errors)

**Everything else is accumulated data**, not a "mode".

### Recommended Approach: Simplified 3-State Model

| State | Detection | Display |
|-------|-----------|---------|
| **Alert** | CVEs > 0 OR errors > 0 OR security issues | Warning + what needs fixing |
| **Swarm** | `activeAgents > 1` | Full swarm details |
| **Normal** | Default | Model + branch + compact data summary |

### Detection Priority
```
Alert → Swarm → Normal
```

### Alternative: Section-Based Display

Instead of modes, show/hide sections based on data existence:
- Always show: Model, branch, context %
- Show if exists: Database stats, cost
- Show if > 0: Patterns, memories
- Show if active: Swarm details
- Show if issues: Alerts

---

## Implementation Notes

### Manual Mode Toggle

No native keyboard shortcut support. Implement via state file:

```bash
# Toggle commands
echo "swarm" > ~/.claude-statusline-mode   # Force swarm view
echo "database" > ~/.claude-statusline-mode # Force database view
echo "auto" > ~/.claude-statusline-mode     # Back to auto-detect
```

In statusline script:
```bash
MODE_OVERRIDE=$(cat ~/.claude-statusline-mode 2>/dev/null)
if [ "$MODE_OVERRIDE" != "auto" ] && [ -n "$MODE_OVERRIDE" ]; then
    MODE="$MODE_OVERRIDE"
else
    # Auto-detect mode
fi
```

### Terminal Width Detection

**Known Issue**: `tput cols` always returns 80 inside Claude Code ([GitHub #5430](https://github.com/anthropics/claude-code/issues/5430) - closed as NOT_PLANNED)

**Working Solution**: Use `stty` via `/dev/tty`:
```bash
# This works! Returns actual terminal width
COLS=$(stty size 2>/dev/null </dev/tty | cut -d' ' -f2 || echo 120)

# Reserve space for Claude Code's right-side messages
USABLE_COLS=$((COLS - 45))
```

| Method | Result | Notes |
|--------|--------|-------|
| `tput cols` | Always 80 | Hardcoded, wrong |
| `stty size </dev/tty` | **Correct** (e.g., 168) | Use this! |
| `$COLUMNS` | Not set | Unavailable |

**Why reserve 45 columns?**

Claude Code displays dynamic messages on the right side:
- `Context left until auto-compact: XX%`
- `Approaching limit`
- IDE integration status

These reduce usable width. The `-45` is a safe buffer (ccstatusline uses `-40`).

**Full-width separator example:**
```bash
COLS=$(stty size 2>/dev/null </dev/tty | cut -d' ' -f2 || echo 120)
printf '━%.0s' $(seq 1 $COLS)
echo
```

**Right-aligned text example:**
```bash
LEFT="Model: Opus"
RIGHT="Context: 42%"
PADDING=$((COLS - ${#LEFT} - ${#RIGHT}))
printf "%s%*s%s\n" "$LEFT" "$PADDING" "" "$RIGHT"
```

### Performance Considerations

- CLI call (`npx @claude-flow/cli@latest hooks statusline --json`) adds ~1-2s latency
- Use timeout to prevent blocking: `timeout 2s npx ...`
- Cache results where possible
- File reads are fast (<10ms)

### Color Codes Reference

```bash
RST="\033[0m"      # Reset
BOLD="\033[1m"     # Bold
DIM="\033[2m"      # Dim
RED="\033[31m"     # Red
GRN="\033[32m"     # Green
YEL="\033[33m"     # Yellow
BLU="\033[34m"     # Blue
MAG="\033[35m"     # Magenta
CYN="\033[36m"     # Cyan
```

### Git History

| Commit | Description |
|--------|-------------|
| `9a35ab7` | Original 311-line comprehensive statusline |
| `e2e9baf` | Added statusline.sh to gitignore |
| `224f264` | Consolidated to single statusline.sh |
| `9b65547` | Adaptive statusline with conditional modes |
| `7222cec` | Hide algorithm line when not active |

---

## References

- [Claude Code Statusline Docs](https://code.claude.com/docs/en/statusline)
- [Claude Code Slash Commands](https://code.claude.com/docs/en/slash-commands)
- [ccstatusline GitHub](https://github.com/sirmalloc/ccstatusline)
- [claude-code-statusline GitHub](https://github.com/rz1989s/claude-code-statusline)
- [Creating The Perfect Claude Code Status Line](https://www.aihero.dev/creating-the-perfect-claude-code-status-line)
