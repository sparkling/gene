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

### From Claude Code stdin (Native)
| Data Point | Source |
|------------|--------|
| Model display name | `model.display_name` |
| Model ID | `model.id` |
| Current directory | `workspace.current_dir` |
| Project directory | `workspace.project_dir` |
| Context used % | `context_window.used_percentage` |
| Context remaining % | `context_window.remaining_percentage` |
| Context window size | `context_window.context_window_size` |
| Total input tokens | `context_window.total_input_tokens` |
| Total output tokens | `context_window.total_output_tokens` |
| Cache creation tokens | `context_window.current_usage.cache_creation_input_tokens` |
| Cache read tokens | `context_window.current_usage.cache_read_input_tokens` |
| Session cost USD | `cost.total_cost_usd` |
| Session duration | `cost.total_duration_ms` |
| API duration | `cost.total_api_duration_ms` |
| Lines added | `cost.total_lines_added` |
| Lines removed | `cost.total_lines_removed` |
| Session ID | `session_id` |
| Claude Code version | `version` |

### From CLI Statusline Hook
| Data Point | Source |
|------------|--------|
| Git branch | `user.gitBranch` |
| Domains completed | `v3Progress.domainsCompleted` |
| Total domains | `v3Progress.totalDomains` |
| DDD progress | `v3Progress.dddProgress` |
| Patterns learned | `v3Progress.patternsLearned` |
| Sessions completed | `v3Progress.sessionsCompleted` |
| Security status | `security.status` |
| CVEs fixed | `security.cvesFixed` |
| Total CVEs | `security.totalCves` |
| Active agents | `swarm.activeAgents` |
| Max agents | `swarm.maxAgents` |
| Coordination active | `swarm.coordinationActive` |
| Memory MB | `system.memoryMB` |
| Context % | `system.contextPct` |
| Intelligence % | `system.intelligencePct` |
| Sub-agents | `system.subAgents` |

### From Intelligence File
| Data Point | Source |
|------------|--------|
| Total patterns | `stats.total_patterns` |
| Total memories | `stats.total_memories` |
| Total trajectories | `stats.total_trajectories` |
| Total errors | `stats.total_errors` |
| Session count | `stats.session_count` |
| Agents count | `agents.length` |
| Edges count | `edges.length` |
| Learning algorithm stats | `learning.stats.*` |

### From Daemon State
| Data Point | Source |
|------------|--------|
| Daemon running | `running` |
| Started at | `startedAt` |
| Worker run count | `workers.*.runCount` |
| Worker success count | `workers.*.successCount` |
| Worker failure count | `workers.*.failureCount` |
| Worker avg duration | `workers.*.averageDurationMs` |
| Worker last run | `workers.*.lastRun` |
| Worker next run | `workers.*.nextRun` |
| Worker is running | `workers.*.isRunning` |

### From Metrics Files
| Data Point | Source |
|------------|--------|
| Flash attention speedup | `performance.json → flashAttention.speedup` |
| HNSW search time | `performance.json → hnsw.searchTimeMs` |
| Routing cost savings | `performance.json → routing.costSavings` |
| Patterns consolidated | `consolidation.json → patternsConsolidated` |
| Memory cleaned | `consolidation.json → memoryCleaned` |

### From Security Audit
| Data Point | Source |
|------------|--------|
| Audit status | `audit-status.json → status` |
| Critical issues | `audit-status.json → issues.critical` |
| High issues | `audit-status.json → issues.high` |
| Medium issues | `audit-status.json → issues.medium` |
| Low issues | `audit-status.json → issues.low` |
| Last scan | `audit-status.json → lastScan` |

### From Git
| Data Point | Command |
|------------|---------|
| Current branch | `git branch --show-current` |
| Uncommitted changes | `git status --porcelain \| wc -l` |
| Last commit | `git log --oneline -1` |
| Stash count | `git stash list \| wc -l` |

### From Database Files
| Data Point | Source |
|------------|--------|
| User DB size | `du -h data/user.db` |
| User vectors (est) | `file_size / 1536` |
| Ops DB size | `du -h data/operational.db` |
| Ops vectors (est) | `file_size / 768` |

**Total: 60+ unique data points**

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
