# Claude Code Configuration - Gene Platform

## Claude-Flow: MCP Only

NEVER use `npx @claude-flow/cli` or Bash for claude-flow operations.

MCP tools handle coordination. Task tool spawns agents. Claude Code tools (Edit, Write) do file work.

**Always run from project root.** Claude-flow creates `.claude-flow/` in the current directory. Running from subdirectories creates duplicate state.

---

## Swarm

**Use when:** 3+ files, new features, refactoring, security, performance

**Skip when:** single file, simple fixes, config, questions

**Pattern:**
1. Spawn ALL agents in ONE message (`run_in_background: true`)
2. Tell user, STOP, WAIT
3. Synthesize when results arrive

**Never:** poll status, add tool calls after spawning

---

## File Organization

Never root. Use: `/src`, `/tests`, `/docs`, `/config`, `/scripts`

---

## Batch Operations

ONE message: TodoWrite, Task tool, file ops, memory ops

---

## Git

Commit after every operation. Don't wait or batch commits.

---

## Reminders

- Do what was asked; nothing more
- Edit existing files over creating new
- Never proactively create docs
