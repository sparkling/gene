---
name: "slt"
description: "Switch between statusline templates. Use /slt list to see options, /slt next or /slt prev to cycle, or /slt {name} to set directly."
invocation: "/slt"
---

# Statusline Switcher

Switch between 24 different statusline display formats organized by purpose.

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

## 24 Available Templates by Category

### Category 1: Minimalist (Focus Mode)
*Purpose: Distraction-free coding, minimal cognitive load*

| Template | Lines | Content |
|----------|-------|---------|
| `zen` | 1 | Model only: `Opus` |
| `focus` | 1 | Model + branch: `Opus │ ⎇ main` |

### Category 2: Developer Daily (Code-Centric)
*Purpose: Active development, git awareness*

| Template | Lines | Content |
|----------|-------|---------|
| `dev` | 1 | `Opus │ gene ⎇ main │ +156/-23 │ 2 uncommitted` |
| `git` | 2 | L1: Model + dir + branch / L2: Changes, stash, last commit |

### Category 3: Knowledge Base (Vector/Learning)
*Purpose: Building RuVector KB, tracking ingestion*

| Template | Lines | Content |
|----------|-------|---------|
| `vectors` | 1 | `Opus │ 📊 819 user │ 682 ops │ 477 memories` |
| `learning` | 2 | L1: Model + patterns / L2: Trajectories, Q-values, sessions |

### Category 4: Operations (Daemon/Workers)
*Purpose: System health monitoring*

| Template | Lines | Content |
|----------|-------|---------|
| `daemon` | 1 | `Opus │ ● daemon │ 383 runs │ 7 workers` |
| `workers` | 2-3 | L1: Model + daemon / L2-3: Per-worker stats |

### Category 5: Cost & Efficiency
*Purpose: Budget awareness, token optimization*

| Template | Lines | Content |
|----------|-------|---------|
| `cost` | 1 | `Opus │ $0.0123 │ 15K in │ 4K out │ 42% ctx` |
| `tokens` | 2 | L1: Cost + context / L2: Cache stats, efficiency |

### Category 6: Security
*Purpose: Security posture awareness*

| Template | Lines | Content |
|----------|-------|---------|
| `secure` | 1 | `Opus │ 🔒 CLEAN │ 0 CVEs │ last scan 2h ago` |
| `audit` | 2 | L1: Status / L2: Issue breakdown by severity |

### Category 7: Swarm (Multi-Agent)
*Purpose: Swarm coordination visibility*

| Template | Lines | Content |
|----------|-------|---------|
| `swarm` | 1 | `Opus │ ⬡ 5/15 agents │ hierarchical-mesh` |
| `agents` | 2-3 | L1: Topology / L2-3: Active agent types |

### Category 8: Performance
*Purpose: Speed and optimization metrics*

| Template | Lines | Content |
|----------|-------|---------|
| `perf` | 1 | `Opus │ ⚡ HNSW 3ms │ Flash 1.0x │ 51.6% saved` |
| `speed` | 2 | L1: Search/attention / L2: Cache, worker durations |

### Category 9: Project Progress
*Purpose: Milestone and domain tracking*

| Template | Lines | Content |
|----------|-------|---------|
| `progress` | 1 | `Opus │ 📈 3/5 domains │ DDD 60% │ 14 sessions` |
| `project` | 2 | L1: Domains / L2: V3 progress bar, patterns |

### Category 10: Dashboard (Comprehensive)
*Purpose: Full visibility, all metrics*

| Template | Lines | Content |
|----------|-------|---------|
| `dashboard` | 3 | L1: Identity / L2: Cost + tokens / L3: Vectors + learning |
| `full` | 5-6 | All categories visible |

### Legacy Templates

| Template | Lines | Content |
|----------|-------|---------|
| `minimal` | 1 | Model, directory, branch |
| `compact` | 1 | Model, branch, context, cost, daemon |
| `adaptive` | 2-4 | Auto-detects swarm/learning/database/idle modes |

## 3-Mode Adaptation System

Every template automatically adapts to these priority modes:

| Mode | Detection | Display Change |
|------|-----------|----------------|
| **Alert** | CVEs > 0 OR critical errors | Red warning header, CVE count |
| **Swarm** | activeAgents > 1 | Green swarm header, agent count |
| **Normal** | Default | Standard category display |

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
.claude/slt.sh dashboard # comprehensive 3-line view
```

## Use Case Quick Reference

| Scenario | Recommended Template |
|----------|---------------------|
| Deep focus work | `zen`, `focus` |
| Daily development | `dev`, `git` |
| Building knowledge base | `vectors`, `learning` |
| Monitoring system | `daemon`, `workers` |
| Watching budget | `cost`, `tokens` |
| Security audit | `secure`, `audit` |
| Multi-agent work | `swarm`, `agents` |
| Optimization work | `perf`, `speed` |
| Project tracking | `progress`, `project` |
| Full overview | `dashboard`, `full` |

## State File

Location: `.claude/statusline-state`

Contains just the template name (no extension):
```
adaptive
```

## Template Directory

Location: `.claude/statuslines/`

Structure:
```
statuslines/
├── lib/
│   └── common.sh    # Shared functions
├── zen.sh           # Minimalist
├── focus.sh         # Minimalist
├── dev.sh           # Developer
├── git.sh           # Developer
├── vectors.sh       # Knowledge Base
├── learning.sh      # Knowledge Base
├── daemon.sh        # Operations
├── workers.sh       # Operations
├── cost.sh          # Cost
├── tokens.sh        # Cost
├── secure.sh        # Security
├── audit.sh         # Security
├── swarm.sh         # Swarm
├── agents.sh        # Swarm
├── perf.sh          # Performance
├── speed.sh         # Performance
├── progress.sh      # Project
├── project.sh       # Project
├── dashboard.sh     # Comprehensive
├── full.sh          # Comprehensive
├── adaptive.sh      # Legacy adaptive
├── compact.sh       # Legacy compact
└── minimal.sh       # Legacy minimal
```

## Notes

- Changes take effect on next statusline render
- No restart required
- Template selection persists across sessions
- Invalid template names show error and available options
- All templates use shared library for consistent behavior
