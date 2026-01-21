#!/bin/bash
#
# Learning Data Migration Script (Simple Version - No SQLite Dependencies)
#
# This script:
# 1. Archives legacy files
# 2. Extracts valid data from intelligence.json to clean JSON
# 3. Generates SQL migration files for later import
# 4. Cleans up garbage data
#
# Usage:
#   ./scripts/migrate-learning-data-simple.sh [--dry-run] [--verbose]
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

# Logging
log() { echo "[MIGRATE] $1"; }
verbose() { $VERBOSE && echo "[DEBUG] $1" || true; }

# Execute command (respects dry-run)
execute() {
  if $DRY_RUN; then
    echo "[DRY-RUN] $*"
  else
    "$@"
  fi
}

log "========================================"
log "Learning Data Migration (Simple)"
log "========================================"
$DRY_RUN && log "DRY RUN MODE - No changes will be made"

# Step 1: Create directories
log ""
log "[Step 1] Creating directories..."
execute mkdir -p "$DATA_DIR"
execute mkdir -p "$ARCHIVE_DIR"

# Step 2: Archive legacy files
log ""
log "[Step 2] Archiving legacy files..."

[ -f "$RUVECTOR_JSON" ] && execute cp "$RUVECTOR_JSON" "$ARCHIVE_DIR/intelligence.json" && log "  Archived intelligence.json"
[ -f "$RUVECTOR_DB" ] && execute cp "$RUVECTOR_DB" "$ARCHIVE_DIR/ruvector.db" && log "  Archived ruvector.db"
[ -f "$SWARM_HNSW" ] && execute cp "$SWARM_HNSW" "$ARCHIVE_DIR/hnsw.index" && log "  Archived hnsw.index"
[ -d "$RUVECTOR_DIR/workers" ] && execute cp -r "$RUVECTOR_DIR/workers" "$ARCHIVE_DIR/workers" && log "  Archived workers/"

# Step 3: Extract and filter data
log ""
log "[Step 3] Extracting valid data..."

if [ -f "$RUVECTOR_JSON" ] && ! $DRY_RUN; then
  node << 'NODEJS_SCRIPT'
const fs = require('fs');
const path = require('path');

const projectRoot = process.cwd();
const dataDir = path.join(projectRoot, 'data');
const ruvectorJson = path.join(projectRoot, '.ruvector/intelligence.json');

try {
  const data = JSON.parse(fs.readFileSync(ruvectorJson, 'utf-8'));

  // Filter valid patterns
  const patterns = Object.entries(data.patterns || {})
    .filter(([k, v]) => v.visits > 0 && v.q_value >= 0 && v.q_value <= 1)
    .map(([k, v]) => ({ key: k, ...v }));

  // Filter valid memories (remove garbage with empty/generic content)
  const seen = new Set();
  const memories = (data.memories || []).filter(m => {
    const hasContent = m.content && m.content.trim().length > 10 && m.content.trim() !== 'Agent:';
    const isUnique = !seen.has(m.id);
    if (isUnique && hasContent) { seen.add(m.id); return true; }
    return false;
  });

  // Filter valid trajectories (deduplicate)
  const trajSeen = new Set();
  const trajectories = (data.trajectories || []).filter(t => {
    const isUnique = !trajSeen.has(t.id);
    if (isUnique && t.state && t.action) { trajSeen.add(t.id); return true; }
    return false;
  });

  // Stats
  const origStats = data.stats || {};
  console.log(`  Original: ${origStats.total_patterns || 0} patterns, ${origStats.total_memories || 0} memories, ${origStats.total_trajectories || 0} trajectories`);
  console.log(`  Filtered: ${patterns.length} patterns, ${memories.length} memories, ${trajectories.length} trajectories`);
  console.log(`  Garbage removed: ${(origStats.total_memories || 0) - memories.length} memories (empty/duplicate)`);

  // Write cleaned data
  const cleanedData = {
    version: '1.0.0',
    migrated_at: new Date().toISOString(),
    source: ruvectorJson,
    stats: {
      patterns: patterns.length,
      memories: memories.length,
      trajectories: trajectories.length
    },
    patterns,
    memories: memories.map(m => ({
      id: m.id,
      memory_type: m.memory_type,
      content: m.content,
      metadata: m.metadata,
      timestamp: m.timestamp
      // Skip embeddings - they're sparse/garbage TF-IDF style
    })),
    trajectories
  };

  fs.writeFileSync(
    path.join(dataDir, 'operational.json'),
    JSON.stringify(cleanedData, null, 2)
  );
  console.log(`  Wrote data/operational.json`);

  // Generate SQL migration script
  const sqlLines = [
    '-- Migration script for operational.db',
    '-- Generated: ' + new Date().toISOString(),
    '',
    '-- Schema',
    'CREATE TABLE IF NOT EXISTS patterns (',
    '  id INTEGER PRIMARY KEY AUTOINCREMENT,',
    '  state TEXT NOT NULL,',
    '  action TEXT NOT NULL,',
    '  q_value REAL DEFAULT 0.5,',
    '  visits INTEGER DEFAULT 0,',
    '  last_update INTEGER,',
    '  created_at INTEGER DEFAULT (strftime(\'%s\', \'now\')),',
    '  UNIQUE(state, action)',
    ');',
    '',
    'CREATE TABLE IF NOT EXISTS memories (',
    '  id TEXT PRIMARY KEY,',
    '  memory_type TEXT NOT NULL,',
    '  content TEXT,',
    '  metadata TEXT,',
    '  timestamp INTEGER,',
    '  created_at INTEGER DEFAULT (strftime(\'%s\', \'now\'))',
    ');',
    '',
    'CREATE TABLE IF NOT EXISTS trajectories (',
    '  id TEXT PRIMARY KEY,',
    '  state TEXT NOT NULL,',
    '  action TEXT NOT NULL,',
    '  outcome TEXT,',
    '  reward REAL,',
    '  timestamp INTEGER,',
    '  created_at INTEGER DEFAULT (strftime(\'%s\', \'now\'))',
    ');',
    '',
    'CREATE INDEX IF NOT EXISTS idx_patterns_state ON patterns(state);',
    'CREATE INDEX IF NOT EXISTS idx_memories_type ON memories(memory_type);',
    'CREATE INDEX IF NOT EXISTS idx_trajectories_state ON trajectories(state);',
    '',
    '-- Data',
    ''
  ];

  // Pattern inserts
  for (const p of patterns) {
    const state = p.state.replace(/'/g, "''");
    const action = p.action.replace(/'/g, "''");
    sqlLines.push(`INSERT OR REPLACE INTO patterns (state, action, q_value, visits, last_update) VALUES ('${state}', '${action}', ${p.q_value}, ${p.visits}, ${p.last_update});`);
  }

  sqlLines.push('');

  // Memory inserts
  for (const m of cleanedData.memories) {
    const content = (m.content || '').replace(/'/g, "''");
    const metadata = JSON.stringify(m.metadata || {}).replace(/'/g, "''");
    sqlLines.push(`INSERT OR REPLACE INTO memories (id, memory_type, content, metadata, timestamp) VALUES ('${m.id}', '${m.memory_type}', '${content}', '${metadata}', ${m.timestamp});`);
  }

  sqlLines.push('');

  // Trajectory inserts
  for (const t of trajectories) {
    const state = t.state.replace(/'/g, "''");
    const action = t.action.replace(/'/g, "''");
    const outcome = (t.outcome || '').replace(/'/g, "''");
    sqlLines.push(`INSERT OR REPLACE INTO trajectories (id, state, action, outcome, reward, timestamp) VALUES ('${t.id}', '${state}', '${action}', '${outcome}', ${t.reward}, ${t.timestamp});`);
  }

  fs.writeFileSync(path.join(dataDir, 'operational.sql'), sqlLines.join('\n'));
  console.log(`  Wrote data/operational.sql`);

  // Write user.sql (empty schema for graph data)
  const userSql = `-- Schema for user.db (graph data)
-- Generated: ${new Date().toISOString()}

CREATE TABLE IF NOT EXISTS nodes (
  id TEXT PRIMARY KEY,
  type TEXT NOT NULL,
  label TEXT,
  properties TEXT,
  created_at INTEGER DEFAULT (strftime('%s', 'now')),
  updated_at INTEGER DEFAULT (strftime('%s', 'now'))
);

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

CREATE INDEX IF NOT EXISTS idx_nodes_type ON nodes(type);
CREATE INDEX IF NOT EXISTS idx_edges_source ON edges(source_id);
CREATE INDEX IF NOT EXISTS idx_edges_target ON edges(target_id);
`;

  fs.writeFileSync(path.join(dataDir, 'user.sql'), userSql);
  console.log(`  Wrote data/user.sql`);

} catch (err) {
  console.error('  Error:', err.message);
  process.exit(1);
}
NODEJS_SCRIPT
fi

# Step 4: Clean up legacy files
log ""
log "[Step 4] Cleaning up legacy files..."

[ -d "$RUVECTOR_DIR" ] && execute rm -rf "$RUVECTOR_DIR" && log "  Removed .ruvector/ (archived)"
[ -f "$RUVECTOR_DB" ] && execute rm -f "$RUVECTOR_DB" && log "  Removed ruvector.db (archived)"

log "  Keeping .swarm/memory.db (actively used by CLI)"
log "  Keeping .swarm/hnsw.index (may be needed for vector search)"

# Summary
log ""
log "========================================"
log "Migration Complete!"
log "========================================"
log ""
log "Output files:"
log "  - $DATA_DIR/operational.json (clean patterns, memories, trajectories)"
log "  - $DATA_DIR/operational.sql (SQL migration script)"
log "  - $DATA_DIR/user.sql (empty schema for graph data)"
log ""
log "Archived files:"
log "  - $ARCHIVE_DIR/"
log ""
log "To create SQLite databases, run:"
log "  sqlite3 data/operational.db < data/operational.sql"
log "  sqlite3 data/user.db < data/user.sql"
log ""
log "Still active (CLI uses these):"
log "  - $SWARM_MEMORY"
log "  - $SWARM_HNSW"

$DRY_RUN && log "" && log "[DRY RUN] No actual changes were made."
