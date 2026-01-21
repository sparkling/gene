#!/bin/bash
#
# Learning Data Migration Script (Bash Version)
#
# Consolidates fragmented learning data from 5 locations into unified structure:
#
# Source (current state):
#   1. .swarm/memory.db (147KB) - SQLite, PRIMARY used by CLI
#   2. .swarm/hnsw.index (1.5MB) - REDB format HNSW index
#   3. .ruvector/intelligence.json (598KB) - Legacy JSON with patterns/memories
#   4. ruvector.db (1.5MB) - REDB format, UNUSED
#   5. data/ - EMPTY placeholder
#
# Target (consolidated):
#   - data/operational.db - Claude Flow internals, patterns, routing
#   - data/user.db - User graph data (reserved for future)
#   - archive/ - Legacy files preserved
#
# Usage:
#   ./scripts/migrate-learning-data.sh [--dry-run] [--verbose]
#

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Parse arguments
DRY_RUN=false
VERBOSE=false
for arg in "$@"; do
  case $arg in
    --dry-run) DRY_RUN=true ;;
    --verbose) VERBOSE=true ;;
  esac
done

# Source paths
SWARM_MEMORY="$PROJECT_ROOT/.swarm/memory.db"
SWARM_HNSW="$PROJECT_ROOT/.swarm/hnsw.index"
RUVECTOR_JSON="$PROJECT_ROOT/.ruvector/intelligence.json"
RUVECTOR_DB="$PROJECT_ROOT/ruvector.db"
RUVECTOR_DIR="$PROJECT_ROOT/.ruvector"

# Target paths
DATA_DIR="$PROJECT_ROOT/data"
ARCHIVE_DIR="$PROJECT_ROOT/archive/migration-$TIMESTAMP"
OPERATIONAL_DB="$DATA_DIR/operational.db"
USER_DB="$DATA_DIR/user.db"

# Logging
log() { echo "[MIGRATE] $1"; }
verbose() { $VERBOSE && echo "[DEBUG] $1" || true; }
warn() { echo "[WARN] $1" >&2; }
error() { echo "[ERROR] $1" >&2; }

# Execute command (respects dry-run)
execute() {
  if $DRY_RUN; then
    echo "[DRY-RUN] $*"
  else
    "$@"
  fi
}

log "========================================"
log "Learning Data Migration Script"
log "========================================"
$DRY_RUN && log "DRY RUN MODE - No changes will be made"

# Step 1: Create directories
log ""
log "[Step 1] Creating target directories..."
execute mkdir -p "$DATA_DIR"
execute mkdir -p "$ARCHIVE_DIR"

# Step 2: Archive legacy files FIRST
log ""
log "[Step 2] Archiving legacy files..."

if [ -f "$RUVECTOR_JSON" ]; then
  log "  Archiving $RUVECTOR_JSON"
  execute cp "$RUVECTOR_JSON" "$ARCHIVE_DIR/intelligence.json"
fi

if [ -f "$RUVECTOR_DB" ]; then
  log "  Archiving $RUVECTOR_DB"
  execute cp "$RUVECTOR_DB" "$ARCHIVE_DIR/ruvector.db"
fi

if [ -f "$SWARM_HNSW" ]; then
  log "  Archiving $SWARM_HNSW"
  execute cp "$SWARM_HNSW" "$ARCHIVE_DIR/hnsw.index"
fi

if [ -d "$RUVECTOR_DIR/workers" ]; then
  log "  Archiving $RUVECTOR_DIR/workers"
  execute cp -r "$RUVECTOR_DIR/workers" "$ARCHIVE_DIR/workers"
fi

# Step 3: Extract and filter valid data from intelligence.json
log ""
log "[Step 3] Extracting valid data from intelligence.json..."

if [ -f "$RUVECTOR_JSON" ]; then
  # Use Node.js to parse JSON and extract valid data
  EXTRACTED_DATA=$(node -e "
    const fs = require('fs');
    const data = JSON.parse(fs.readFileSync('$RUVECTOR_JSON', 'utf-8'));

    // Extract valid patterns
    const patterns = Object.entries(data.patterns || {})
      .filter(([k, v]) => v.visits > 0 && v.q_value >= 0 && v.q_value <= 1)
      .map(([k, v]) => ({ key: k, ...v }));

    // Extract valid memories (filter garbage with empty content)
    const seen = new Set();
    const memories = (data.memories || [])
      .filter(m => {
        const hasContent = m.content && m.content.trim().length > 10;
        const isUnique = !seen.has(m.id);
        if (isUnique && hasContent) { seen.add(m.id); return true; }
        return false;
      });

    // Extract valid trajectories (deduplicate)
    const trajSeen = new Set();
    const trajectories = (data.trajectories || [])
      .filter(t => {
        const isUnique = !trajSeen.has(t.id);
        if (isUnique && t.state && t.action) { trajSeen.add(t.id); return true; }
        return false;
      });

    console.log(JSON.stringify({
      stats: {
        original: data.stats || {},
        filtered: {
          patterns: patterns.length,
          memories: memories.length,
          trajectories: trajectories.length
        }
      },
      patterns,
      memories,
      trajectories
    }));
  " 2>/dev/null || echo '{"error":"parse_failed"}')

  if echo "$EXTRACTED_DATA" | grep -q '"error"'; then
    warn "Failed to parse intelligence.json"
  else
    ORIG_PATTERNS=$(echo "$EXTRACTED_DATA" | node -e "const d=require('fs').readFileSync(0,'utf-8');console.log(JSON.parse(d).stats.original.total_patterns||0)")
    ORIG_MEMORIES=$(echo "$EXTRACTED_DATA" | node -e "const d=require('fs').readFileSync(0,'utf-8');console.log(JSON.parse(d).stats.original.total_memories||0)")
    ORIG_TRAJECTORIES=$(echo "$EXTRACTED_DATA" | node -e "const d=require('fs').readFileSync(0,'utf-8');console.log(JSON.parse(d).stats.original.total_trajectories||0)")

    FILT_PATTERNS=$(echo "$EXTRACTED_DATA" | node -e "const d=require('fs').readFileSync(0,'utf-8');console.log(JSON.parse(d).stats.filtered.patterns)")
    FILT_MEMORIES=$(echo "$EXTRACTED_DATA" | node -e "const d=require('fs').readFileSync(0,'utf-8');console.log(JSON.parse(d).stats.filtered.memories)")
    FILT_TRAJECTORIES=$(echo "$EXTRACTED_DATA" | node -e "const d=require('fs').readFileSync(0,'utf-8');console.log(JSON.parse(d).stats.filtered.trajectories)")

    log "  Original: $ORIG_PATTERNS patterns, $ORIG_MEMORIES memories, $ORIG_TRAJECTORIES trajectories"
    log "  After filtering: $FILT_PATTERNS patterns, $FILT_MEMORIES memories, $FILT_TRAJECTORIES trajectories"
    log "  Garbage removed: $((ORIG_MEMORIES - FILT_MEMORIES)) memories"

    # Save filtered data for import
    if ! $DRY_RUN; then
      echo "$EXTRACTED_DATA" > "$DATA_DIR/extracted_data.json"
    fi
  fi
else
  warn "intelligence.json not found"
fi

# Step 4: Create operational database schema
log ""
log "[Step 4] Creating operational database..."

SCHEMA_SQL="
-- Patterns table: Q-learning state-action values
CREATE TABLE IF NOT EXISTS patterns (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  state TEXT NOT NULL,
  action TEXT NOT NULL,
  q_value REAL DEFAULT 0.5,
  visits INTEGER DEFAULT 0,
  last_update INTEGER,
  created_at INTEGER DEFAULT (strftime('%s', 'now')),
  UNIQUE(state, action)
);

-- Memories table: Semantic memories with embeddings
CREATE TABLE IF NOT EXISTS memories (
  id TEXT PRIMARY KEY,
  memory_type TEXT NOT NULL,
  content TEXT,
  embedding BLOB,
  metadata TEXT,
  timestamp INTEGER,
  created_at INTEGER DEFAULT (strftime('%s', 'now'))
);

-- Trajectories table: Learning trajectories for SONA
CREATE TABLE IF NOT EXISTS trajectories (
  id TEXT PRIMARY KEY,
  state TEXT NOT NULL,
  action TEXT NOT NULL,
  outcome TEXT,
  reward REAL,
  timestamp INTEGER,
  created_at INTEGER DEFAULT (strftime('%s', 'now'))
);

-- HNSW index metadata
CREATE TABLE IF NOT EXISTS hnsw_metadata (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  namespace TEXT NOT NULL,
  dimension INTEGER NOT NULL,
  m INTEGER DEFAULT 16,
  ef_construction INTEGER DEFAULT 200,
  total_vectors INTEGER DEFAULT 0,
  last_rebuild INTEGER,
  created_at INTEGER DEFAULT (strftime('%s', 'now'))
);

-- Sessions table
CREATE TABLE IF NOT EXISTS sessions (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  session_id TEXT UNIQUE,
  started_at INTEGER,
  ended_at INTEGER,
  patterns_learned INTEGER DEFAULT 0,
  memories_stored INTEGER DEFAULT 0,
  trajectories_recorded INTEGER DEFAULT 0
);

-- Migration metadata
CREATE TABLE IF NOT EXISTS _migrations (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  migration_name TEXT NOT NULL,
  executed_at INTEGER DEFAULT (strftime('%s', 'now')),
  source_files TEXT,
  records_migrated INTEGER
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_patterns_state ON patterns(state);
CREATE INDEX IF NOT EXISTS idx_memories_type ON memories(memory_type);
CREATE INDEX IF NOT EXISTS idx_trajectories_state ON trajectories(state);
CREATE INDEX IF NOT EXISTS idx_trajectories_timestamp ON trajectories(timestamp);
"

if ! $DRY_RUN; then
  # Use Node.js with sql.js (WASM SQLite) since sqlite3 CLI is not available
  node -e "
    const initSqlJs = require('sql.js');
    const fs = require('fs');

    (async () => {
      const SQL = await initSqlJs();
      const db = new SQL.Database();

      // Create schema
      db.run(\`$SCHEMA_SQL\`);

      // Load extracted data if available
      const dataPath = '$DATA_DIR/extracted_data.json';
      if (fs.existsSync(dataPath)) {
        const data = JSON.parse(fs.readFileSync(dataPath, 'utf-8'));

        // Insert patterns
        const patternStmt = db.prepare('INSERT OR REPLACE INTO patterns (state, action, q_value, visits, last_update) VALUES (?, ?, ?, ?, ?)');
        for (const p of data.patterns || []) {
          patternStmt.run([p.state, p.action, p.q_value, p.visits, p.last_update]);
        }
        patternStmt.free();

        // Insert memories (skip embeddings for now - they're sparse/garbage)
        const memStmt = db.prepare('INSERT OR REPLACE INTO memories (id, memory_type, content, metadata, timestamp) VALUES (?, ?, ?, ?, ?)');
        for (const m of data.memories || []) {
          memStmt.run([m.id, m.memory_type, m.content, JSON.stringify(m.metadata), m.timestamp]);
        }
        memStmt.free();

        // Insert trajectories
        const trajStmt = db.prepare('INSERT OR REPLACE INTO trajectories (id, state, action, outcome, reward, timestamp) VALUES (?, ?, ?, ?, ?, ?)');
        for (const t of data.trajectories || []) {
          trajStmt.run([t.id, t.state, t.action, t.outcome, t.reward, t.timestamp]);
        }
        trajStmt.free();

        // Record migration
        db.run('INSERT INTO _migrations (migration_name, source_files, records_migrated) VALUES (?, ?, ?)',
          ['initial_consolidation', '$RUVECTOR_JSON', (data.patterns?.length || 0) + (data.memories?.length || 0) + (data.trajectories?.length || 0)]);

        console.log('Data imported successfully');
      }

      // Save database
      const dbData = db.export();
      const buffer = Buffer.from(dbData);
      fs.writeFileSync('$OPERATIONAL_DB', buffer);
      db.close();

      console.log('Created $OPERATIONAL_DB');
    })();
  " 2>&1 || {
    warn "sql.js not available, creating empty database marker"
    touch "$OPERATIONAL_DB.pending"
  }
fi

# Step 5: Create user database schema
log ""
log "[Step 5] Creating user database..."

USER_SCHEMA_SQL="
-- User graph nodes
CREATE TABLE IF NOT EXISTS nodes (
  id TEXT PRIMARY KEY,
  type TEXT NOT NULL,
  label TEXT,
  properties TEXT,
  created_at INTEGER DEFAULT (strftime('%s', 'now')),
  updated_at INTEGER DEFAULT (strftime('%s', 'now'))
);

-- User graph edges
CREATE TABLE IF NOT EXISTS edges (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  source_id TEXT NOT NULL,
  target_id TEXT NOT NULL,
  relationship TEXT NOT NULL,
  weight REAL DEFAULT 1.0,
  properties TEXT,
  created_at INTEGER DEFAULT (strftime('%s', 'now')),
  FOREIGN KEY (source_id) REFERENCES nodes(id),
  FOREIGN KEY (target_id) REFERENCES nodes(id)
);

-- User embeddings
CREATE TABLE IF NOT EXISTS user_embeddings (
  id TEXT PRIMARY KEY,
  node_id TEXT,
  embedding BLOB,
  dimension INTEGER,
  created_at INTEGER DEFAULT (strftime('%s', 'now')),
  FOREIGN KEY (node_id) REFERENCES nodes(id)
);

-- Indexes
CREATE INDEX IF NOT EXISTS idx_nodes_type ON nodes(type);
CREATE INDEX IF NOT EXISTS idx_edges_source ON edges(source_id);
CREATE INDEX IF NOT EXISTS idx_edges_target ON edges(target_id);
CREATE INDEX IF NOT EXISTS idx_edges_relationship ON edges(relationship);
"

if ! $DRY_RUN; then
  node -e "
    const initSqlJs = require('sql.js');
    const fs = require('fs');

    (async () => {
      const SQL = await initSqlJs();
      const db = new SQL.Database();
      db.run(\`$USER_SCHEMA_SQL\`);
      const dbData = db.export();
      fs.writeFileSync('$USER_DB', Buffer.from(dbData));
      db.close();
      console.log('Created $USER_DB');
    })();
  " 2>&1 || {
    warn "sql.js not available, creating empty database marker"
    touch "$USER_DB.pending"
  }
fi

# Step 6: Clean up legacy files
log ""
log "[Step 6] Cleaning up legacy files..."

if [ -d "$RUVECTOR_DIR" ]; then
  log "  Removing $RUVECTOR_DIR (archived)"
  execute rm -rf "$RUVECTOR_DIR"
fi

if [ -f "$RUVECTOR_DB" ]; then
  log "  Removing $RUVECTOR_DB (archived)"
  execute rm -f "$RUVECTOR_DB"
fi

log "  Keeping .swarm/memory.db (actively used by CLI)"
log "  Keeping .swarm/hnsw.index (may be needed for vector search)"

# Step 7: Clean up temp files
log ""
log "[Step 7] Cleaning up temporary files..."
if [ -f "$DATA_DIR/extracted_data.json" ]; then
  execute rm -f "$DATA_DIR/extracted_data.json"
fi

# Summary
log ""
log "========================================"
log "Migration Complete!"
log "========================================"
log ""
log "New database locations:"
log "  - $OPERATIONAL_DB (patterns, memories, trajectories)"
log "  - $USER_DB (user graph data - empty, ready for use)"
log ""
log "Archived files:"
log "  - $ARCHIVE_DIR/"
log ""
log "Still active (CLI uses these):"
log "  - $SWARM_MEMORY"
log "  - $SWARM_HNSW"

if $DRY_RUN; then
  log ""
  log "[DRY RUN] No actual changes were made. Run without --dry-run to execute."
fi
