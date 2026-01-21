# Architecture Decision Records

**Document ID:** ADR-INDEX
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

This directory contains Architecture Decision Records (ADRs) documenting significant technical decisions for the Gene platform. ADRs follow the MADR 3.0 format and are integrated with Claude Flow's ReasoningBank for pattern learning.

---

## ADR Index

| ADR | Title | Status | Date |
|-----|-------|--------|------|
| [ADR-001](./ADR-001-three-database-architecture.md) | Three-Database Architecture | Accepted | Jan 2026 |
| [ADR-002](./ADR-002-knowledge-graph-engine.md) | Knowledge Graph Engine Selection | Accepted | Jan 2026 |
| [ADR-003](./ADR-003-embedding-model-selection.md) | Embedding Model Selection | Accepted | Jan 2026 |
| [ADR-004](./ADR-004-federation-strategy.md) | SPARQL Federation Strategy | Accepted | Jan 2026 |
| [ADR-005](./ADR-005-claude-flow-memory.md) | Claude Flow Memory Configuration | Accepted | Jan 2026 |
| [ADR-006](./ADR-006-hardware-optimized-hnsw.md) | Hardware-Optimized HNSW Parameters | Accepted | Jan 2026 |
| [ADR-007](./ADR-007-automated-learning-configuration.md) | Automated Learning Configuration | Accepted | Jan 2026 |
| [ADR-008](./ADR-008-mcp-server-management.md) | MCP Server Management | Accepted | Jan 2026 |

---

## Decision Categories

| Category | ADRs | Description |
|----------|------|-------------|
| **Architecture** | ADR-001 | System structure and database topology |
| **Data** | ADR-002, ADR-004 | Data storage, retrieval, and federation |
| **AI/ML** | ADR-003, ADR-005, ADR-007 | Embeddings, vectors, and learning systems |
| **Infrastructure** | ADR-006, ADR-008 | Hardware optimization and MCP server management |

---

## ADR Format (MADR 3.0)

All ADRs follow this structure:

```markdown
# ADR-XXX: {Title}

## Status
{Proposed | Accepted | Deprecated | Superseded by ADR-XXX}

## Context
What is the issue motivating this decision?

## Decision
What is the change we're proposing/doing?

## Consequences
### Positive
### Negative
### Neutral

## Options Considered
### Option 1: {Name}
### Option 2: {Name}

## Related Decisions
## References
```

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [44-ARCHITECTURE](../40-product/44-ARCHITECTURE.md) | Implements these decisions |
| [43-DATA-SOURCES](../40-product/data-sources/43-00-INDEX.md) | Data source inventory |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Initial ADR index with 5 decisions |
| 1.1 | January 2026 | Engineering | Added ADR-006, ADR-007, ADR-008 for hardware, learning, MCP |
