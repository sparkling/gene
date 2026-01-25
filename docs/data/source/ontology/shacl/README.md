---
id: ontology-shacl
title: "SHACL Shapes"
type: ontology
parent: ../_index.md
last_updated: 2026-01-25
status: active
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

## Files

| File | Description |
|------|-------------|
| [datasource-shapes.ttl](./datasource-shapes.ttl) | Complete validation shapes for data sources |

## Namespace

```turtle
@prefix shapes: <https://gene.example.org/shapes/datasource#> .
```

## Shapes

| Shape | Target Class | Required Properties |
|-------|--------------|---------------------|
| DataSourceShape | ds:DataSource | id, title, category, subcategory, status |
| CategoryShape | ds:Category | id, title, classification |
| SubcategoryShape | ds:Subcategory | id, title, category |
| VersionShape | ds:Version | versionNumber, releaseDate |
| LicenseShape | ds:License | licenseName |
| AccessMethodShape | ds:AccessMethod | accessMethodType |
| SchemaShape | ds:Schema | (fields recommended) |
| FieldShape | ds:Field | fieldName, fieldType |
| IdentifierShape | ds:Identifier | identifierName |
| CrossReferenceShape | ds:CrossReference | referencesSource |

## Severity Levels

| Level | Severity | Description |
|-------|----------|-------------|
| Critical | `sh:Violation` | Must fix - data is invalid |
| Warning | `sh:Warning` | Should fix - recommended |
| Info | `sh:Info` | Nice to have - optional |

## Key Constraints

### DataSourceShape

| Property | Constraint | Severity |
|----------|------------|----------|
| id | 2-50 chars, pattern `^[a-z0-9][a-z0-9.-]*[a-z0-9]$` | Violation |
| title | 2-100 chars | Violation |
| status | one of: Draft, Active, Deprecated | Violation |
| tier | one of: Tier1, Tier2, Tier3 | Warning |
| description | min 50 chars | Warning |
| hasLicense | required for active sources | Violation |

### VersionShape

| Property | Constraint | Severity |
|----------|------------|----------|
| versionNumber | required, string | Violation |
| releaseDate | required, xsd:date | Violation |
| checksum | pattern `^(sha256|sha512|md5):[a-f0-9]+$` | Info |

### LicenseShape

| Property | Constraint | Severity |
|----------|------------|----------|
| licenseName | required, min 2 chars | Violation |
| spdxIdentifier | pattern `^[A-Za-z0-9.-]+$` | Info |
| commercialUseAllowed | boolean | Warning |

## Validation

```bash
# Using pySHACL
pyshacl -s datasource-shapes.ttl -e ../owl/datasource-ontology.ttl data.ttl

# Using Apache Jena
shacl validate --shapes datasource-shapes.ttl --data data.ttl
```

## SPARQL Constraints

The shapes include SPARQL-based constraints:

- Active sources must have license information
- Non-draft sources should have at least one version
- Subcategories should have at least one data source
