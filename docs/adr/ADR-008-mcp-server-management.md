# ADR-008: MCP Server Management

## Status
Accepted

## Context

Claude Flow runs as an MCP (Model Context Protocol) server, providing tools for:
- Agent spawning and coordination
- Memory storage and retrieval
- Swarm orchestration
- Learning and pattern recognition

Managing the MCP server requires understanding:
1. How to restart the server (to apply configuration changes)
2. How to set runtime properties
3. How to configure environment variables
4. When each method is appropriate

## Decision

### Restart Methods

| Method | Command | Resets Uptime | Reloads Config | Preserves Session |
|--------|---------|---------------|----------------|-------------------|
| **System Reset** | `system_reset({ component: "all", confirm: true })` | ✅ Yes | ✅ Yes | ✅ Yes |
| MCP Reconnect | `/mcp reconnect claude-flow` | ❌ No | ❌ No | ✅ Yes |
| Full Restart | Exit + relaunch Claude Code | ✅ Yes | ✅ Yes | ❌ No |
| Daemon Restart | `daemon stop` + `daemon start` | N/A | Partial | ✅ Yes |

**Primary Restart Method**: Use `system_reset` for most cases.

```javascript
// Recommended restart method
mcp__claude-flow__system_reset({
  component: "all",  // "all", "metrics", "agents", "tasks"
  confirm: true
})
```

### Property Setting Methods

| Method | When to Use | Persistence |
|--------|-------------|-------------|
| `config_set` MCP tool | Runtime changes | Until reset |
| `.claude-flow/config.json` edit | Permanent changes | Always |
| `.mcp.json` env vars | Startup config (heap, threads) | Always |

**Primary Method**: Use `config_set` for immediate changes.

```javascript
// Set runtime config
mcp__claude-flow__config_set({ key: "hnsw.m", value: 48 })

// Get current value
mcp__claude-flow__config_get({ key: "hnsw.m" })

// List all config
mcp__claude-flow__config_list({})

// Import multiple settings
mcp__claude-flow__config_import({
  config: { "hnsw.m": 48, "learning.autoLearn": true },
  merge: true
})
```

### Environment Variables (.mcp.json)

For settings that must be set at startup (heap size, worker threads):

```json
{
  "mcpServers": {
    "claude-flow": {
      "command": "npx",
      "args": ["claude-flow@v3alpha", "mcp", "start"],
      "env": {
        "CLAUDE_FLOW_MAX_AGENTS": "100",
        "CLAUDE_FLOW_WORKER_THREADS": "16",
        "NODE_OPTIONS": "--max-old-space-size=65536"
      }
    }
  }
}
```

**Note**: Changes to `.mcp.json` require full Claude Code restart.

### Verification

```javascript
// Check uptime (should be near 0 after restart)
mcp__claude-flow__system_status({ verbose: true })

// Verify config applied
mcp__claude-flow__config_get({ key: "hnsw.m" })

// Check learning systems active
mcp__claude-flow__hooks_intelligence({ showStatus: true })
```

## Consequences

### Positive
- Clear hierarchy of restart methods
- `system_reset` preserves session while reloading config
- MCP tools provide immediate feedback
- Environment variables for startup-only settings

### Negative
- `.mcp.json` changes require full restart
- Multiple configuration sources can conflict
- Need to know which method for which change

### Neutral
- Runtime config in `.claude-flow/config.json`
- Startup config in `.mcp.json`

## Configuration Sources (Priority Order)

1. **MCP tool calls** (immediate, session-scoped)
2. **`.claude-flow/config.json`** (persistent, reloaded on system_reset)
3. **`.mcp.json` env vars** (startup only, requires full restart)

## Common Operations

### Apply Config Changes
```javascript
mcp__claude-flow__config_set({ key: "...", value: ... })
mcp__claude-flow__system_reset({ component: "all", confirm: true })
```

### Change Heap Size
1. Edit `.mcp.json`: `NODE_OPTIONS": "--max-old-space-size=65536"`
2. Exit Claude Code
3. Relaunch Claude Code

### Debug Config Issues
```javascript
mcp__claude-flow__config_export({})  // See all settings
mcp__claude-flow__system_info({})    // System information
mcp__claude-flow__system_health({ deep: true })  // Health check
```

## Related Decisions
- [ADR-006](./ADR-006-hardware-optimized-hnsw.md): Hardware-Optimized HNSW Parameters
- [ADR-007](./ADR-007-automated-learning-configuration.md): Automated Learning Configuration

## References
- [MCP-CONFIGURATION-GUIDE.md](../architecture/MCP-CONFIGURATION-GUIDE.md)
- [Model Context Protocol Specification](https://modelcontextprotocol.io/)
