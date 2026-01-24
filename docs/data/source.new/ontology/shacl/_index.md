---
id: ontology-shacl
title: "SHACL Shapes"
type: ontology
parent: ../_index.md
last_updated: 2026-01-23
status: placeholder
tags: [shacl, validation, constraints, data-quality]
---

# SHACL Shapes

**Parent:** [Ontology](../_index.md)

SHACL (Shapes Constraint Language) definitions for data validation, quality constraints, and structural requirements.

## Purpose

- Validate data source metadata against required schema
- Enforce cardinality and datatype constraints
- Define required vs optional properties
- Enable automated data quality checking

## Planned Files

| File | Description |
|------|-------------|
| `data-source-shapes.ttl` | Validation shapes for data source entries |
| `schema-shapes.ttl` | Shapes for schema documentation |
| `relationship-shapes.ttl` | Constraints on entity relationships |

## Example Constraints

```turtle
# Example: Data source must have title and category
ex:DataSourceShape a sh:NodeShape ;
    sh:targetClass :DataSource ;
    sh:property [
        sh:path :title ;
        sh:minCount 1 ;
        sh:datatype xsd:string ;
    ] ;
    sh:property [
        sh:path :category ;
        sh:minCount 1 ;
        sh:in ( :genetics :compounds :diseases :pathways ) ;
    ] .
```

## Validation Levels

| Level | Severity | Description |
|-------|----------|-------------|
| Critical | sh:Violation | Must fix before import |
| Warning | sh:Warning | Should fix, non-blocking |
| Info | sh:Info | Recommendation only |
