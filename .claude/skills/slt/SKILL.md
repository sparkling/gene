---
name: "slt"
description: "Switch between statusline templates. Use /slt list to see options, /slt next or /slt prev to cycle, or /slt {name} to set directly."
invocation: "/slt"
---

# Statusline Switcher

Switch between 8 curated statusline templates organized by purpose.

## Commands

**List templates:**
```
/slt
/slt list
```

**Cycle templates:**
```
/slt next     # Next template in list
/slt prev     # Previous template in list
/slt n        # Shortcut for next
/slt p        # Shortcut for prev
```

**Set specific template:**
```
/slt <template-name>
```

## 8 Templates by Category

### Focus (Minimal)
*Purpose: Distraction-free coding, minimal cognitive load*

| Template | Lines | Content |
|----------|-------|---------|
| `zen` | 1 | Model only with mode indicators: `Opus` or `⬡ Opus (3 agents)` |

### Development
*Purpose: Active development, git awareness*

| Template | Lines | Content |
|----------|-------|---------|
| `dev` | 1 | `Opus │ gene ⎇ main │ +156/-23 │ 2 uncommitted` |

### Cost & Tokens
*Purpose: Budget awareness, token optimization*

| Template | Lines | Content |
|----------|-------|---------|
| `cost` | 1 | `Opus │ $0.0123 │ 15K in │ 4K out │ 42% ctx` |

### Operations
*Purpose: System health monitoring*

| Template | Lines | Content |
|----------|-------|---------|
| `daemon` | 1 | `Opus │ ● daemon │ 383 runs │ 7 workers` |

### Security
*Purpose: Security posture awareness*

| Template | Lines | Content |
|----------|-------|---------|
| `secure` | 1 | `Opus │ 🔒 CLEAN │ 0 CVEs │ last scan 2h ago` |

### Coordination
*Purpose: Multi-agent swarm visibility*

| Template | Lines | Content |
|----------|-------|---------|
| `swarm` | 1 | `Opus │ ⬡ 3 running │ hierarchical │ main +2` |

### Dashboard
*Purpose: Full visibility, most actionable metrics*

| Template | Lines | Content |
|----------|-------|---------|
| `full` | 3 | L1: Model + git / L2: Context + cost + tokens / L3: Agents + daemon |

### Adaptive
*Purpose: Auto-selects display based on current activity*

| Template | Lines | Content |
|----------|-------|---------|
| `adaptive` | 2-3 | Auto-detects alert/swarm/active/idle modes |

## 3-Mode Adaptation System

Every template automatically adapts to these priority modes:

| Mode | Detection | Display Change |
|------|-----------|----------------|
| **Alert** | CVEs > 0 OR critical errors | Red warning header, CVE count |
| **Swarm** | Task tool agents > 1 | Green swarm header, running agent count |
| **Normal** | Default | Standard category display |

## Data Sources (All Verified Accurate)

| Data Point | Source | Update Trigger |
|------------|--------|----------------|
| Model name | stdin (Claude Code) | Per render |
| Git branch/changes | `git status`, `git diff` | Per render |
| Cost/tokens | stdin (Claude Code) | Per render |
| Context % | stdin (Claude Code) | Per render |
| Running agents | `/tmp/claude/.../tasks/*.output` | Per render |
| Registry agents | `.claude-flow/agents/store.json` | MCP operations |
| Daemon status | `.claude-flow/daemon-state.json` | Daemon lifecycle |
| Security/CVEs | `.claude-flow/metrics/audit-status.json` | Security scan |
| Topology | `claude-flow.config.json` | Config change |

## Removed Data (Previously Unreliable)

The following were removed as they showed fake/static data:

- Vector counts (were calculated from file sizes, not actual counts)
- V3 progress percentages (internal dev metrics with schema mismatches)
- HNSW latency (no real data source)
- Flash Attention speedup (no real measurement)
- Learning algorithm stats (sparse/missing data)
- Pattern/trajectory counts (often 0 or stale)

## Quick Execution

```bash
.claude/slt.sh [ARG]
```

Where ARG is: `list`, `n`, `p`, `next`, `prev`, or a template name.

**Examples:**
```bash
.claude/slt.sh          # list templates
.claude/slt.sh n        # next
.claude/slt.sh p        # prev
.claude/slt.sh zen      # minimal focus mode
.claude/slt.sh full     # comprehensive 3-line view
```

## Use Case Quick Reference

| Scenario | Recommended Template |
|----------|---------------------|
| Deep focus work | `zen` |
| Daily development | `dev` |
| Watching budget | `cost` |
| Monitoring system | `daemon` |
| Security audit | `secure` |
| Multi-agent work | `swarm` |
| Full overview | `full` |
| Auto-detect needs | `adaptive` |

## State File

Location: `.claude/statusline-state`

Contains just the template name (no extension):
```
full
```

## Template Directory

Location: `.claude/statuslines/`

Structure:
```
statuslines/
├── lib/
│   └── common.sh    # Shared functions (data loading, mode detection)
├── zen.sh           # Focus (minimal)
├── dev.sh           # Development
├── cost.sh          # Cost & Tokens
├── daemon.sh        # Operations
├── secure.sh        # Security
├── swarm.sh         # Coordination
├── full.sh          # Dashboard (comprehensive)
└── adaptive.sh      # Adaptive (context-aware)
```

## Agent Count Terminology

| Term | Meaning | Source |
|------|---------|--------|
| **running** | Real Task tool agents executing now | `/tmp/claude/.../tasks/` |
| **registered** | MCP registry entries (may not be active) | `.claude-flow/agents/store.json` |
| **queued** | Tasks waiting for agent assignment | `.claude-flow/tasks.json` |

## Notes

- Changes take effect on next statusline render
- No restart required
- Template selection persists across sessions
- Invalid template names show error and available options
- All templates use shared library for consistent behavior
- All displayed data is verified to update dynamically
