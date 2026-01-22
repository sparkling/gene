# Statusline Data Accuracy Report

**Date**: 2026-01-21
**Auditor**: Code Analyzer Agent
**Scope**: `/home/claude/src/gene/.claude/statuslines/` data sources

---

## Executive Summary

| Data Source | Accuracy | Freshness | Status |
|-------------|----------|-----------|--------|
| `agents/store.json` | HIGH | CURRENT | Working correctly |
| `daemon-state.json` | MEDIUM | CURRENT | Field name mismatch found |
| `v3-progress.json` | LOW | STALE | Schema mismatch, outdated |
| `intelligence.json` | MEDIUM | STALE | Data is ~38 hours old |

---

## 1. `.claude-flow/agents/store.json`

### File Status
- **Last Modified**: 2026-01-21 16:18:49 (recent)
- **Size**: ~3KB
- **Update Frequency**: Dynamic (updated on agent spawn/terminate)

### Content Analysis
```json
{
  "agents": { /* 8 agents registered */ },
  "version": "3.0.0"
}
```

**Actual Data**:
- 8 agents registered
- Agent types: `coder` (3), `researcher` (2), `reviewer` (1), undefined (2)
- All agents status: `idle`
- All agents health: `1` (healthy)
- All taskCount: `0`

### lib/common.sh Reading (load_swarm function)

| JQ Query | Expected | Actual | CORRECT |
|----------|----------|--------|---------|
| `.agents \| length` | 8 | 8 | YES |
| `.agents \| to_entries \| .[].value.agentType` | agent types | coder, researcher, coder, researcher, coder, reviewer | YES |
| `[.agents \| to_entries \| .[].value.taskCount] \| add` | 0 | 0 | YES |

### Issues Found
- **NONE** - This data source is being read correctly

---

## 2. `.claude-flow/daemon-state.json`

### File Status
- **Last Modified**: 2026-01-21 16:25:44 (very recent - 7 minutes ago)
- **Size**: ~3KB
- **Update Frequency**: Dynamic (updated every worker run)

### Content Analysis
```json
{
  "running": true,
  "startedAt": "2026-01-20T09:33:44.611Z",
  "workers": {
    "map": { "runCount": 328, "successCount": 328, "isRunning": false },
    "audit": { "runCount": 333, "successCount": 0, "failureCount": 333, "isRunning": false },
    "optimize": { "runCount": 245, "successCount": 0, "failureCount": 245, "isRunning": false },
    "consolidate": { "runCount": 163, "successCount": 163, "isRunning": false },
    "testgaps": { "runCount": 198, "successCount": 1, "isRunning": true }
  },
  "savedAt": "2026-01-21T15:25:44.941Z"
}
```

### lib/common.sh Reading (load_daemon function)

| JQ Query | Expected | Actual | CORRECT |
|----------|----------|--------|---------|
| `.running` | true | true | YES |
| `[.workers[].runCount] \| add` | 1267 | 1267 | YES |
| `[.workers[].successCount] \| add` | 492 | (null - no output) | **NO** |
| `[.workers[] \| select(.running == true)] \| length` | 1 | 0 | **NO** |

### Issues Found

**CRITICAL BUG**: Field name mismatch in `load_daemon`:

1. **`.running` vs `.isRunning`**: The daemon-state.json uses `isRunning` but the jq query looks for `.running`:
   ```bash
   # Line 129 in common.sh:
   WORKERS_RUNNING=$(jq_safe "$DAEMON_FILE" '[.workers[] | select(.running == true)] | length' "0")
   ```
   **Actual field name in JSON**: `isRunning` (not `running`)

   **Result**: `testgaps` worker shows as running (`isRunning: true`) but the statusline reports 0 running workers.

2. **Null handling in successCount**: Some workers have `successCount: 0` which when combined, jq may return null instead of integer if any worker lacks the field.

### Recommended Fix
Change line 129 in common.sh:
```bash
# BEFORE:
WORKERS_RUNNING=$(jq_safe "$DAEMON_FILE" '[.workers[] | select(.running == true)] | length' "0")

# AFTER:
WORKERS_RUNNING=$(jq_safe "$DAEMON_FILE" '[.workers[] | select(.isRunning == true)] | length' "0")
```

---

## 3. `.claude-flow/metrics/v3-progress.json`

### File Status
- **Last Modified**: 2026-01-20 03:11:09 (~37 hours old)
- **Size**: ~250 bytes
- **Update Frequency**: Manual/CLI (`progress sync`)

### Content Analysis
```json
{
  "domains": { "completed": 3, "total": 5 },
  "ddd": { "progress": 60 },
  "swarm": { "activeAgents": 0, "maxAgents": 50 },
  "memory": { "entries": 12, "sizeMb": 18 },
  "patternsLearned": 140,
  "sessionsCompleted": 14
}
```

### lib/common.sh Reading (load_progress function)

| JQ Query | Expected | Actual | CORRECT |
|----------|----------|--------|---------|
| `.overall.percentage` | v3 progress % | **null** (field missing) | **NO** |
| `.domains.completed` | 3 | 3 | YES |
| `.domains.total` | 5 | 5 | YES |
| `.ddd.percentage` | DDD progress | **null** (`.ddd.progress` exists) | **NO** |

### Issues Found

**SCHEMA MISMATCH**: The common.sh expects different field names than what's in the file:

1. **Missing `.overall.percentage`**: The file has no `overall` key at all. The script expects:
   ```bash
   V3_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.overall.percentage // 0' "0")
   ```
   But the file structure has no `overall.percentage`.

2. **Field name `.ddd.percentage` vs `.ddd.progress`**:
   ```bash
   DDD_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.ddd.percentage // 0' "0")
   ```
   Actual field is `.ddd.progress` (value: 60), not `.ddd.percentage`.

3. **Stale Data**: File is 37+ hours old, not being updated dynamically.

4. **Swarm data mismatch**: Shows `activeAgents: 0` but agents/store.json shows 8 agents registered.

### Recommended Fixes
1. Update common.sh to match actual schema:
   ```bash
   # Fix DDD progress field name
   DDD_PROGRESS=$(jq_safe "$PROGRESS_FILE" '.ddd.progress // 0' "0")
   ```
2. Either add `overall.percentage` to v3-progress.json or calculate from domains:
   ```bash
   # Calculate from domains
   V3_PROGRESS=$(jq_safe "$PROGRESS_FILE" '(.domains.completed / .domains.total * 100) | floor' "0")
   ```

---

## 4. `.ruvector/intelligence.json`

### File Status
- **Last Modified**: 2026-01-20 02:19:46 (~38 hours old)
- **Size**: 584KB (large file with embeddings)
- **Update Frequency**: Sporadic (last pattern update: 1768871846 = 2026-01-20 02:17:26)

### Content Analysis
```json
{
  "patterns": { /* 2 patterns */ },
  "memories": [ /* 477 memories with embeddings */ ],
  "trajectories": [ /* 293 trajectories */ ],
  "stats": {
    "total_patterns": 2,
    "total_memories": 477,
    "total_trajectories": 293,
    "total_errors": 0,
    "session_count": 7,
    "last_session": 1768871732
  }
}
```

### lib/common.sh Reading (load_intel function)

| JQ Query | Expected | Actual | CORRECT |
|----------|----------|--------|---------|
| `.stats.total_patterns` | 2 | 2 | YES |
| `.stats.total_memories` | 477 | 477 | YES |
| `.stats.total_trajectories` | 293 | 293 | YES |
| `.stats.total_errors` | 0 | 0 | YES |
| `.stats.session_count` | 7 | 7 | YES |
| `.learning.stats...` | algorithm data | **null** (no learning key) | **NO** |

### Issues Found

1. **Missing `.learning` key**: The `load_intel` function tries to read best algorithm data:
   ```bash
   RV_BEST_ALGO=$(jq -r '.learning.stats // {} | to_entries | ...')
   ```
   But the file has no `learning` key. Result: `RV_BEST_ALGO="none"`, `RV_CONVERGENCE=0`, `RV_Q_AVG=0`.

2. **Stale data**: The file hasn't been updated in ~38 hours. The `last_update` timestamps in patterns are from 1768871846 (2026-01-20 02:17:26).

3. **Pattern Q-values exist but in different location**: The patterns have `q_value` directly on them, but common.sh looks for `.learning.stats`.

### Recommended Fix
Either:
1. Update common.sh to read Q-values from patterns directly:
   ```bash
   RV_BEST_ALGO=$(jq -r '.patterns | to_entries | sort_by(-.value.q_value) | .[0].key // "none"' "$INTEL_FILE" 2>/dev/null || echo "none")
   ```
2. Or ensure the intelligence system writes a `.learning.stats` section.

---

## 5. Missing Data Sources

### `audit-status.json`
- **Status**: FILE DOES NOT EXIST
- **Expected Path**: `.claude-flow/metrics/audit-status.json`
- **Impact**: Security statuslines will always show "unknown" status, 0 CVEs
- **Used By**: `detect_mode()`, `load_security()`, `secure.sh`, `audit.sh`

### `performance.json`
- **Status**: EXISTS but minimal data
- **Path**: `.claude-flow/metrics/performance.json`
- **Content**:
  ```json
  {
    "flashAttention": { "enabled": true, "speedup": "1.0x" },
    "hnsw": { "searchTimeMs": 3, "indexSize": 12 },
    "routing": { "costSavings": "51.6%" }
  }
  ```
- **Issues**:
  - `speedup` has format "1.0x" - common.sh strips trailing 'x' (handled correctly)
  - Missing `.cache.hitRate` (defaults to 0)
  - Missing `.memory.savedPercent` (defaults to 0)

### `claude-flow.config.json`
- **Status**: FILE DOES NOT EXIST
- **Expected Path**: `./claude-flow.config.json`
- **Impact**: Swarm topology/max agents will use defaults ("hierarchical", 15)
- **Used By**: `load_swarm()`

### Vector Databases (`data/*.db`)
- **Status**: EMPTY DIRECTORY
- **Path**: `./data/`
- **Content**: Only `.gitkeep` file exists
- **Impact**: `load_vectors()` will report 0 vectors, 0K size for all DBs

---

## Summary of Required Fixes

### Critical (Breaking Functionality)

| File | Issue | Fix Required |
|------|-------|--------------|
| `lib/common.sh:129` | `.running` should be `.isRunning` | Change field name |
| `lib/common.sh:278` | `.ddd.percentage` should be `.ddd.progress` | Change field name |

### Medium (Incorrect/Missing Data)

| Issue | Impact | Recommendation |
|-------|--------|----------------|
| Missing `.overall.percentage` in v3-progress.json | V3_PROGRESS always 0 | Calculate from domains or add field |
| Missing `.learning.stats` in intelligence.json | Best algo always "none" | Read from patterns directly or update writer |
| Stale v3-progress.json (37h old) | Incorrect progress displayed | Run `progress sync` regularly |
| Stale intelligence.json (38h old) | Learning data outdated | Trigger learning updates on activity |

### Low (Enhancement Opportunities)

| Issue | Recommendation |
|-------|----------------|
| Missing audit-status.json | Create via `security scan` command |
| Missing claude-flow.config.json | Create via `init` command |
| Empty data/ directory | Vector DBs created on first use |

---

## Test Commands

To verify fixes, run these commands:

```bash
# Test daemon worker count (should match isRunning)
jq '[.workers[] | select(.isRunning == true)] | length' .claude-flow/daemon-state.json

# Test v3 progress calculation
jq '(.domains.completed / .domains.total * 100) | floor' .claude-flow/metrics/v3-progress.json

# Test DDD progress with correct field
jq '.ddd.progress' .claude-flow/metrics/v3-progress.json

# Test best algorithm from patterns
jq '.patterns | to_entries | sort_by(-.value.q_value) | .[0].key' .ruvector/intelligence.json
```

---

## Appendix: Data Freshness Timeline

| Data Source | Last Modified | Age | Should Be |
|-------------|---------------|-----|-----------|
| agents/store.json | 2026-01-21 16:18 | ~9 min | Dynamic - OK |
| daemon-state.json | 2026-01-21 16:25 | ~2 min | Dynamic - OK |
| v3-progress.json | 2026-01-20 03:11 | ~37h | Should update more frequently |
| intelligence.json | 2026-01-20 02:19 | ~38h | Should update on learning activity |
