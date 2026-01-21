# MCP Server Configuration Guide

> How to configure, restart, and manage the Claude Flow MCP server

**Last Updated:** 2026-01-22

---

## Table of Contents

1. [Quick Reference](#quick-reference)
2. [How to Restart MCP Server](#how-to-restart-mcp-server)
3. [How to Set MCP Properties](#how-to-set-mcp-properties)
4. [MCP Configuration File (.mcp.json)](#mcp-configuration-file)
5. [Environment Variables](#environment-variables)
6. [Verification](#verification)
7. [Troubleshooting](#troubleshooting)

---

## Quick Reference

| Action | Command |
|--------|---------|
| Restart MCP (preserve session) | `mcp__claude-flow__system_reset({ component: "all", confirm: true })` |
| Set config property | `mcp__claude-flow__config_set({ key: "...", value: ... })` |
| Get config property | `mcp__claude-flow__config_get({ key: "..." })` |
| List all config | `mcp__claude-flow__config_list({})` |
| Check uptime | `mcp__claude-flow__system_status({ verbose: true })` |
| Reconnect MCP | `/mcp reconnect claude-flow` |

---

## How to Restart MCP Server

### Method 1: System Reset (Recommended)

Resets internal state, reloads all configurations, resets uptime to 0. **Does NOT lose MCP connection.**

```javascript
// Via MCP tool
mcp__claude-flow__system_reset({
  component: "all",    // Options: "all", "metrics", "agents", "tasks"
  confirm: true        // Required confirmation
})
```

**Result:** Uptime resets to 0h 0m, all configs reloaded, session preserved.

### Method 2: MCP Reconnect

Reconnects to MCP server without restarting it. **Does NOT reset uptime or reload configs.**

```
/mcp reconnect claude-flow
```

**Use case:** Fix connection issues, not for config changes.

### Method 3: Full Restart (Session Lost)

For changes to `.mcp.json` environment variables:

1. Exit Claude Code completely (`Ctrl+C` or `/exit`)
2. Relaunch Claude Code
3. MCP servers reinitialize with new environment variables

**Use case:** After modifying `.mcp.json` env vars.

### Method 4: Daemon Restart

Restart background workers only (not the MCP server itself):

```bash
# Via CLI (if needed)
npx @claude-flow/cli@latest daemon stop
npx @claude-flow/cli@latest daemon start
```

### Comparison Table

| Method | Resets Uptime | Reloads Config | Preserves Session | Use Case |
|--------|---------------|----------------|-------------------|----------|
| system_reset | ✅ Yes | ✅ Yes | ✅ Yes | **Best for most cases** |
| /mcp reconnect | ❌ No | ❌ No | ✅ Yes | Connection issues |
| Full restart | ✅ Yes | ✅ Yes | ❌ No | `.mcp.json` env changes |
| Daemon restart | N/A | Partial | ✅ Yes | Worker issues |

---

## How to Set MCP Properties

### Method 1: MCP Tools (Recommended)

Use `config_set` for immediate application:

```javascript
// Set boolean
mcp__claude-flow__config_set({
  key: "learning.autoLearn",
  value: true
})

// Set numeric
mcp__claude-flow__config_set({
  key: "hnsw.m",
  value: 48
})

// Set string
mcp__claude-flow__config_set({
  key: "swarm.topology",
  value: "hierarchical-mesh"
})

// Set nested value
mcp__claude-flow__config_set({
  key: "daemon.workers.ultralearn.enabled",
  value: true
})
```

### Method 2: Get Current Values

```javascript
// Get single value
mcp__claude-flow__config_get({ key: "hnsw.m" })

// List all values
mcp__claude-flow__config_list({})

// List with prefix filter
mcp__claude-flow__config_list({ prefix: "hnsw" })

// Export full config
mcp__claude-flow__config_export({})
```

### Method 3: Import Configuration

```javascript
// Import multiple settings at once
mcp__claude-flow__config_import({
  config: {
    "hnsw.m": 48,
    "hnsw.efConstruction": 400,
    "learning.autoLearn": true
  },
  merge: true  // true = merge with existing, false = replace all
})
```

### Method 4: Reset to Defaults

```javascript
// Reset specific key
mcp__claude-flow__config_reset({ key: "hnsw.m" })

// Reset all (caution!)
mcp__claude-flow__config_reset({})
```

---

## MCP Configuration File

The `.mcp.json` file in your project root configures MCP servers.

### Location

```
/home/claude/src/gene/.mcp.json
```

### Structure

```json
{
  "mcpServers": {
    "claude-flow": {
      "command": "npx",
      "args": ["claude-flow@v3alpha", "mcp", "start"],
      "env": {
        "CLAUDE_FLOW_MAX_AGENTS": "100",
        "CLAUDE_FLOW_HNSW_M": "48",
        "CLAUDE_FLOW_FLASH_ATTENTION": "true",
        "NODE_OPTIONS": "--max-old-space-size=65536"
      }
    }
  }
}
```

### Fields

| Field | Description |
|-------|-------------|
| `command` | Executable to run (usually `npx`) |
| `args` | Arguments to pass (package and subcommand) |
| `env` | Environment variables passed to the MCP server |

### Applying .mcp.json Changes

Changes to `.mcp.json` require a **full restart** of Claude Code:

1. Save `.mcp.json`
2. Exit Claude Code
3. Relaunch Claude Code

The `system_reset` MCP tool does NOT reload `.mcp.json` - it only reloads runtime config from `.claude-flow/config.json`.

---

## Environment Variables

### In .mcp.json

Set environment variables in the `env` section:

```json
{
  "mcpServers": {
    "claude-flow": {
      "env": {
        "CLAUDE_FLOW_MAX_AGENTS": "100",
        "CLAUDE_FLOW_HNSW_M": "48",
        "CLAUDE_FLOW_HNSW_EF_CONSTRUCTION": "400",
        "CLAUDE_FLOW_HNSW_EF_SEARCH": "200",
        "CLAUDE_FLOW_FLASH_ATTENTION": "true",
        "CLAUDE_FLOW_WORKER_THREADS": "16",
        "CLAUDE_FLOW_MEMORY_MAX_ENTRIES": "50000000",
        "CLAUDE_FLOW_EMBEDDINGS_CACHE_SIZE": "8192",
        "CLAUDE_FLOW_DAEMON_MAX_CONCURRENT": "16",
        "CLAUDE_FLOW_NEURAL_ENABLED": "true",
        "CLAUDE_FLOW_LORA_RANK": "16",
        "NODE_OPTIONS": "--max-old-space-size=65536"
      }
    }
  }
}
```

### Common Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `CLAUDE_FLOW_MAX_AGENTS` | Maximum concurrent agents | `100` |
| `CLAUDE_FLOW_HNSW_M` | HNSW connections per node | `48` |
| `CLAUDE_FLOW_HNSW_EF_CONSTRUCTION` | HNSW build quality | `400` |
| `CLAUDE_FLOW_HNSW_EF_SEARCH` | HNSW search quality | `200` |
| `CLAUDE_FLOW_FLASH_ATTENTION` | Enable Flash Attention | `true` |
| `CLAUDE_FLOW_WORKER_THREADS` | Worker thread count | `16` |
| `CLAUDE_FLOW_MEMORY_MAX_ENTRIES` | Max memory entries | `50000000` |
| `CLAUDE_FLOW_NEURAL_ENABLED` | Enable neural learning | `true` |
| `NODE_OPTIONS` | Node.js options | `--max-old-space-size=65536` |

### Environment Variables vs Runtime Config

| Aspect | Environment Variables (.mcp.json) | Runtime Config (config_set) |
|--------|-----------------------------------|----------------------------|
| When applied | On MCP server start | Immediately |
| Requires restart | Full Claude Code restart | system_reset only |
| Persistence | Always (file-based) | Until reset |
| Use case | Initial setup, heap size | Dynamic tuning |

---

## Verification

### Check MCP Server Status

```javascript
// Get system status including uptime
mcp__claude-flow__system_status({ verbose: true })
```

Example output:
```
{
  "status": "operational",
  "uptime": "0h 0m 15s",
  "components": {
    "memory": "healthy",
    "neural": "active",
    "swarm": "ready"
  }
}
```

### Verify Configuration Applied

```javascript
// Check specific value
mcp__claude-flow__config_get({ key: "hnsw.m" })

// List all with prefix
mcp__claude-flow__config_list({ prefix: "hnsw" })
```

### Check Intelligence Systems

```javascript
mcp__claude-flow__hooks_intelligence({ showStatus: true })
```

---

## Troubleshooting

### Config Changes Not Taking Effect

1. **For runtime config (config_set):**
   ```javascript
   mcp__claude-flow__system_reset({ component: "all", confirm: true })
   ```

2. **For .mcp.json environment variables:**
   - Exit Claude Code completely
   - Relaunch Claude Code

### MCP Connection Issues

```
/mcp reconnect claude-flow
```

### Check What's Currently Set

```javascript
// Export current config
mcp__claude-flow__config_export({})

// Check system info
mcp__claude-flow__system_info({})
```

### Uptime Not Resetting

Only `system_reset` with `component: "all"` resets uptime:

```javascript
// This DOES reset uptime
mcp__claude-flow__system_reset({ component: "all", confirm: true })

// These do NOT reset uptime
mcp__claude-flow__swarm_shutdown({ graceful: true })  // ❌
/mcp reconnect claude-flow                             // ❌
```

### View MCP Server Logs

```javascript
// Get recent metrics
mcp__claude-flow__system_metrics({ timeRange: "1h" })

// Check health
mcp__claude-flow__system_health({ deep: true })
```

---

## Summary

| Task | Best Method |
|------|-------------|
| Change runtime settings | `config_set` + `system_reset` |
| Change heap size/threads | Edit `.mcp.json` + full restart |
| Fix connection issues | `/mcp reconnect` |
| Verify restart worked | Check uptime via `system_status` |
| Debug config issues | `config_export` + `config_list` |

---

*See also:*
- [RUVECTOR-CONFIGURATION-GUIDE.md](./RUVECTOR-CONFIGURATION-GUIDE.md) - Hardware optimization
- [AUTOMATED-LEARNING-GUIDE.md](./AUTOMATED-LEARNING-GUIDE.md) - Learning settings
