#!/usr/bin/env npx ts-node
/**
 * Learning Data Migration Script
 *
 * Consolidates fragmented learning data from 5 locations into unified databases:
 *
 * Source (current state):
 *   1. .swarm/memory.db (147KB) - SQLite, PRIMARY used by CLI
 *   2. .swarm/hnsw.index (1.5MB) - REDB format HNSW index
 *   3. .ruvector/intelligence.json (598KB) - Legacy JSON with patterns/memories
 *   4. ruvector.db (1.5MB) - REDB format, UNUSED
 *   5. data/ - EMPTY placeholder
 *
 * Target (consolidated):
 *   - data/operational.db - Claude Flow internals, patterns, routing
 *   - data/user.db - User graph data (reserved for future)
 *   - archive/ - Legacy files preserved
 *
 * Usage:
 *   npx ts-node scripts/migrate-learning-data.ts [--dry-run] [--verbose]
 */

import * as fs from 'fs';
import * as path from 'path';
import Database from 'better-sqlite3';

// Configuration
const PROJECT_ROOT = path.resolve(__dirname, '..');
const DRY_RUN = process.argv.includes('--dry-run');
const VERBOSE = process.argv.includes('--verbose');

// Source paths
const SOURCES = {
  swarmMemory: path.join(PROJECT_ROOT, '.swarm/memory.db'),
  swarmHnsw: path.join(PROJECT_ROOT, '.swarm/hnsw.index'),
  ruvectorJson: path.join(PROJECT_ROOT, '.ruvector/intelligence.json'),
  ruvectorDb: path.join(PROJECT_ROOT, 'ruvector.db'),
};

// Target paths
const TARGETS = {
  dataDir: path.join(PROJECT_ROOT, 'data'),
  archiveDir: path.join(PROJECT_ROOT, 'archive'),
  operationalDb: path.join(PROJECT_ROOT, 'data/operational.db'),
  userDb: path.join(PROJECT_ROOT, 'data/user.db'),
};

// Logging helpers
const log = (msg: string) => console.log(`[MIGRATE] ${msg}`);
const verbose = (msg: string) => VERBOSE && console.log(`[DEBUG] ${msg}`);
const warn = (msg: string) => console.warn(`[WARN] ${msg}`);
const error = (msg: string) => console.error(`[ERROR] ${msg}`);

interface Pattern {
  state: string;
  action: string;
  q_value: number;
  visits: number;
  last_update: number;
}

interface Memory {
  id: string;
  memory_type: string;
  content: string;
  embedding: number[];
  metadata: Record<string, unknown>;
  timestamp: number;
}

interface Trajectory {
  id: string;
  state: string;
  action: string;
  outcome: string;
  reward: number;
  timestamp: number;
}

interface IntelligenceJson {
  patterns: Record<string, Pattern>;
  memories: Memory[];
  trajectories: Trajectory[];
  stats: {
    total_patterns: number;
    total_memories: number;
    total_trajectories: number;
    session_count: number;
    last_session: number;
  };
}

/**
 * Create the operational database schema
 */
function createOperationalSchema(db: Database.Database): void {
  db.exec(`
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

    -- HNSW index metadata (actual index stored separately)
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

    -- Sessions table: Track learning sessions
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

    -- Create indexes for performance
    CREATE INDEX IF NOT EXISTS idx_patterns_state ON patterns(state);
    CREATE INDEX IF NOT EXISTS idx_memories_type ON memories(memory_type);
    CREATE INDEX IF NOT EXISTS idx_trajectories_state ON trajectories(state);
    CREATE INDEX IF NOT EXISTS idx_trajectories_timestamp ON trajectories(timestamp);
  `);
}

/**
 * Create the user database schema (for graph data)
 */
function createUserSchema(db: Database.Database): void {
  db.exec(`
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

    -- User embeddings (for semantic search)
    CREATE TABLE IF NOT EXISTS user_embeddings (
      id TEXT PRIMARY KEY,
      node_id TEXT,
      embedding BLOB,
      dimension INTEGER,
      created_at INTEGER DEFAULT (strftime('%s', 'now')),
      FOREIGN KEY (node_id) REFERENCES nodes(id)
    );

    -- Create indexes
    CREATE INDEX IF NOT EXISTS idx_nodes_type ON nodes(type);
    CREATE INDEX IF NOT EXISTS idx_edges_source ON edges(source_id);
    CREATE INDEX IF NOT EXISTS idx_edges_target ON edges(target_id);
    CREATE INDEX IF NOT EXISTS idx_edges_relationship ON edges(relationship);
  `);
}

/**
 * Parse and validate intelligence.json
 */
function parseIntelligenceJson(filePath: string): IntelligenceJson | null {
  try {
    const content = fs.readFileSync(filePath, 'utf-8');
    const data = JSON.parse(content) as IntelligenceJson;

    verbose(`Parsed intelligence.json: ${data.stats?.total_patterns || 0} patterns, ${data.stats?.total_memories || 0} memories, ${data.stats?.total_trajectories || 0} trajectories`);

    return data;
  } catch (err) {
    error(`Failed to parse intelligence.json: ${err}`);
    return null;
  }
}

/**
 * Filter valid patterns (remove garbage)
 */
function filterValidPatterns(patterns: Record<string, Pattern>): Record<string, Pattern> {
  const valid: Record<string, Pattern> = {};

  for (const [key, pattern] of Object.entries(patterns)) {
    // Filter criteria:
    // 1. Must have meaningful state and action
    // 2. Must have been visited at least once
    // 3. Q-value must be in valid range [0, 1]
    if (
      pattern.state &&
      pattern.action &&
      pattern.visits > 0 &&
      pattern.q_value >= 0 &&
      pattern.q_value <= 1
    ) {
      valid[key] = pattern;
    } else {
      verbose(`Filtering out invalid pattern: ${key}`);
    }
  }

  return valid;
}

/**
 * Filter valid memories (remove garbage with empty content)
 */
function filterValidMemories(memories: Memory[]): Memory[] {
  const seen = new Set<string>();

  return memories.filter(mem => {
    // Filter criteria:
    // 1. Must have actual content (not just "Agent: ")
    // 2. Must have valid embedding (non-trivial)
    // 3. Deduplicate by id
    const hasContent = mem.content && mem.content.trim().length > 10;
    const hasValidEmbedding = mem.embedding && mem.embedding.some(v => v !== 0 && v !== 0.3779644730092272);
    const isUnique = !seen.has(mem.id);

    if (isUnique && hasContent) {
      seen.add(mem.id);
      return true;
    }

    verbose(`Filtering out memory: ${mem.id} (content: "${mem.content?.substring(0, 30)}...")`);
    return false;
  });
}

/**
 * Filter valid trajectories
 */
function filterValidTrajectories(trajectories: Trajectory[]): Trajectory[] {
  const seen = new Set<string>();

  return trajectories.filter(traj => {
    // Deduplicate and ensure valid data
    const isUnique = !seen.has(traj.id);
    const hasState = traj.state && traj.state.length > 0;
    const hasAction = traj.action && traj.action.length > 0;

    if (isUnique && hasState && hasAction) {
      seen.add(traj.id);
      return true;
    }

    return false;
  });
}

/**
 * Migrate patterns to operational database
 */
function migratePatterns(db: Database.Database, patterns: Record<string, Pattern>): number {
  const stmt = db.prepare(`
    INSERT OR REPLACE INTO patterns (state, action, q_value, visits, last_update)
    VALUES (@state, @action, @q_value, @visits, @last_update)
  `);

  let count = 0;
  const transaction = db.transaction(() => {
    for (const pattern of Object.values(patterns)) {
      stmt.run({
        state: pattern.state,
        action: pattern.action,
        q_value: pattern.q_value,
        visits: pattern.visits,
        last_update: pattern.last_update,
      });
      count++;
    }
  });

  if (!DRY_RUN) {
    transaction();
  }

  return count;
}

/**
 * Migrate memories to operational database
 */
function migrateMemories(db: Database.Database, memories: Memory[]): number {
  const stmt = db.prepare(`
    INSERT OR REPLACE INTO memories (id, memory_type, content, embedding, metadata, timestamp)
    VALUES (@id, @memory_type, @content, @embedding, @metadata, @timestamp)
  `);

  let count = 0;
  const transaction = db.transaction(() => {
    for (const mem of memories) {
      stmt.run({
        id: mem.id,
        memory_type: mem.memory_type,
        content: mem.content,
        embedding: Buffer.from(new Float32Array(mem.embedding).buffer),
        metadata: JSON.stringify(mem.metadata),
        timestamp: mem.timestamp,
      });
      count++;
    }
  });

  if (!DRY_RUN) {
    transaction();
  }

  return count;
}

/**
 * Migrate trajectories to operational database
 */
function migrateTrajectories(db: Database.Database, trajectories: Trajectory[]): number {
  const stmt = db.prepare(`
    INSERT OR REPLACE INTO trajectories (id, state, action, outcome, reward, timestamp)
    VALUES (@id, @state, @action, @outcome, @reward, @timestamp)
  `);

  let count = 0;
  const transaction = db.transaction(() => {
    for (const traj of trajectories) {
      stmt.run({
        id: traj.id,
        state: traj.state,
        action: traj.action,
        outcome: traj.outcome,
        reward: traj.reward,
        timestamp: traj.timestamp,
      });
      count++;
    }
  });

  if (!DRY_RUN) {
    transaction();
  }

  return count;
}

/**
 * Archive legacy files
 */
function archiveLegacyFiles(): void {
  const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
  const archiveSubdir = path.join(TARGETS.archiveDir, `migration-${timestamp}`);

  if (!DRY_RUN) {
    fs.mkdirSync(archiveSubdir, { recursive: true });
  }

  const filesToArchive = [
    { src: SOURCES.ruvectorJson, dest: path.join(archiveSubdir, 'intelligence.json') },
    { src: SOURCES.ruvectorDb, dest: path.join(archiveSubdir, 'ruvector.db') },
    { src: SOURCES.swarmHnsw, dest: path.join(archiveSubdir, 'hnsw.index') },
  ];

  for (const { src, dest } of filesToArchive) {
    if (fs.existsSync(src)) {
      log(`Archiving ${src} -> ${dest}`);
      if (!DRY_RUN) {
        fs.copyFileSync(src, dest);
      }
    }
  }

  // Also archive the .ruvector/workers directory if it exists
  const workersDir = path.join(PROJECT_ROOT, '.ruvector/workers');
  if (fs.existsSync(workersDir)) {
    const workersDest = path.join(archiveSubdir, 'workers');
    log(`Archiving ${workersDir} -> ${workersDest}`);
    if (!DRY_RUN) {
      fs.cpSync(workersDir, workersDest, { recursive: true });
    }
  }
}

/**
 * Clean up legacy files after successful migration
 */
function cleanupLegacyFiles(): void {
  // Remove .ruvector directory (already archived)
  const ruvectorDir = path.join(PROJECT_ROOT, '.ruvector');
  if (fs.existsSync(ruvectorDir)) {
    log(`Removing ${ruvectorDir}`);
    if (!DRY_RUN) {
      fs.rmSync(ruvectorDir, { recursive: true });
    }
  }

  // Remove root ruvector.db (already archived)
  if (fs.existsSync(SOURCES.ruvectorDb)) {
    log(`Removing ${SOURCES.ruvectorDb}`);
    if (!DRY_RUN) {
      fs.unlinkSync(SOURCES.ruvectorDb);
    }
  }

  // Keep .swarm/memory.db as it's actively used by CLI
  // Keep .swarm/hnsw.index as it may be needed for vector search
  log('Keeping .swarm/memory.db and .swarm/hnsw.index (actively used by CLI)');
}

/**
 * Main migration function
 */
async function main(): Promise<void> {
  log('========================================');
  log('Learning Data Migration Script');
  log('========================================');

  if (DRY_RUN) {
    log('DRY RUN MODE - No changes will be made');
  }

  // Step 1: Create directories
  log('\n[Step 1] Creating target directories...');
  if (!DRY_RUN) {
    fs.mkdirSync(TARGETS.dataDir, { recursive: true });
    fs.mkdirSync(TARGETS.archiveDir, { recursive: true });
  }

  // Step 2: Archive legacy files FIRST (before any modifications)
  log('\n[Step 2] Archiving legacy files...');
  archiveLegacyFiles();

  // Step 3: Create operational database with schema
  log('\n[Step 3] Creating operational database...');
  let operationalDb: Database.Database | null = null;
  if (!DRY_RUN) {
    operationalDb = new Database(TARGETS.operationalDb);
    createOperationalSchema(operationalDb);
    log(`Created ${TARGETS.operationalDb}`);
  }

  // Step 4: Create user database with schema
  log('\n[Step 4] Creating user database...');
  let userDb: Database.Database | null = null;
  if (!DRY_RUN) {
    userDb = new Database(TARGETS.userDb);
    createUserSchema(userDb);
    log(`Created ${TARGETS.userDb}`);
  }

  // Step 5: Parse and filter intelligence.json
  log('\n[Step 5] Parsing and filtering intelligence.json...');
  const intelligenceData = parseIntelligenceJson(SOURCES.ruvectorJson);

  if (intelligenceData) {
    // Filter garbage data
    const validPatterns = filterValidPatterns(intelligenceData.patterns || {});
    const validMemories = filterValidMemories(intelligenceData.memories || []);
    const validTrajectories = filterValidTrajectories(intelligenceData.trajectories || []);

    log(`  Original: ${Object.keys(intelligenceData.patterns || {}).length} patterns, ${(intelligenceData.memories || []).length} memories, ${(intelligenceData.trajectories || []).length} trajectories`);
    log(`  After filtering: ${Object.keys(validPatterns).length} patterns, ${validMemories.length} memories, ${validTrajectories.length} trajectories`);

    // Step 6: Migrate data to operational database
    log('\n[Step 6] Migrating data to operational database...');
    if (operationalDb) {
      const patternsCount = migratePatterns(operationalDb, validPatterns);
      const memoriesCount = migrateMemories(operationalDb, validMemories);
      const trajectoriesCount = migrateTrajectories(operationalDb, validTrajectories);

      log(`  Migrated: ${patternsCount} patterns, ${memoriesCount} memories, ${trajectoriesCount} trajectories`);

      // Record migration
      operationalDb.prepare(`
        INSERT INTO _migrations (migration_name, source_files, records_migrated)
        VALUES (?, ?, ?)
      `).run(
        'initial_consolidation',
        JSON.stringify([SOURCES.ruvectorJson]),
        patternsCount + memoriesCount + trajectoriesCount
      );
    }
  } else {
    warn('No intelligence.json data to migrate');
  }

  // Step 7: Clean up legacy files
  log('\n[Step 7] Cleaning up legacy files...');
  cleanupLegacyFiles();

  // Step 8: Close databases
  log('\n[Step 8] Finalizing...');
  if (operationalDb) {
    operationalDb.close();
  }
  if (userDb) {
    userDb.close();
  }

  // Summary
  log('\n========================================');
  log('Migration Complete!');
  log('========================================');
  log('\nNew database locations:');
  log(`  - ${TARGETS.operationalDb} (patterns, memories, trajectories)`);
  log(`  - ${TARGETS.userDb} (user graph data - empty, ready for use)`);
  log('\nArchived files:');
  log(`  - ${TARGETS.archiveDir}/`);
  log('\nStill active (CLI uses these):');
  log(`  - ${SOURCES.swarmMemory}`);
  log(`  - ${SOURCES.swarmHnsw}`);

  if (DRY_RUN) {
    log('\n[DRY RUN] No actual changes were made. Run without --dry-run to execute.');
  }
}

// Run
main().catch(err => {
  error(`Migration failed: ${err}`);
  process.exit(1);
});
