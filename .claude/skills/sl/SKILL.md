---
name: "Statusline Switcher"
description: "Switch between statusline templates. Use /sl list to see options, /sl next or /sl prev to cycle, or /sl {name} to set directly. Available templates: minimal, compact, adaptive, full."
---

# Statusline Switcher

Switch between different statusline display formats.

## Available Templates

| Template | Lines | Focus | Use Case |
|----------|-------|-------|----------|
| `minimal` | 1 | Identity | Model, directory, branch only |
| `compact` | 1 | Overview | Key metrics on one line |
| `adaptive` | 2-4 | Auto-detect | Mode-based (swarm/learning/database/idle) |
| `full` | 5-6 | Everything | All categories always visible |

## Commands

**List templates:**
```
/sl
/sl list
```

**Cycle templates:**
```
/sl next     # Next template in list
/sl prev     # Previous template in list
/sl n        # Shortcut for next
/sl p        # Shortcut for prev
```

**Set specific template:**
```
/sl minimal
/sl compact
/sl adaptive
/sl full
```

## Quick Execution

**Run this single command:**
```bash
.claude/sl.sh [ARG]
```

Where ARG is: `list`, `n`, `p`, `next`, `prev`, or a template name.

**Examples:**
```bash
.claude/sl.sh          # list templates
.claude/sl.sh n        # next
.claude/sl.sh p        # prev
.claude/sl.sh minimal  # set specific
```

## State File

Location: `.claude/statusline-state`

Contains just the template name (no extension):
```
adaptive
```

## Template Directory

Location: `.claude/statuslines/`

Files:
- `adaptive.sh` - Auto-detects swarm/learning/database/idle modes
- `compact.sh` - Single line with model, branch, daemon, cost
- `full.sh` - Multi-line comprehensive display
- `minimal.sh` - Just model and git branch

## Example Session

```
User: /sl
Claude: Current statusline: adaptive

Available templates:
  minimal   - Identity only (1 line)
  compact   - Key metrics (1 line)
* adaptive  - Mode-based (2-4 lines)
  full      - Everything (5-6 lines)

User: /sl compact
Claude: Statusline switched to: compact

User: /sl next
Claude: Statusline switched to: full

User: /sl prev
Claude: Statusline switched to: compact
```

## Notes

- Changes take effect on next statusline render
- No restart required
- Template selection persists across sessions
- Invalid template names show error and available options
