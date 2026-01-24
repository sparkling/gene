# Global Schema Registry

This directory contains the centralized schema definitions and validation artifacts for the GENE resource collection. These files define the shared data structures, field specifications, and validation rules used across all resource domains.

## Contents

### [shared-fields-registry.json](shared-fields-registry.json)
The authoritative registry of all shared field definitions used across resource schemas. Contains:
- Field names and their data types
- Validation rules and constraints
- Field descriptions and usage guidelines
- Cross-reference mappings between domains

### [data-dictionary.json](data-dictionary.json)
Machine-readable data dictionary containing structured metadata for all fields:
- Complete field specifications
- Type definitions and enumerations
- Relationship mappings
- Required vs optional field designations

### [data-dictionary.md](data-dictionary.md)
Human-readable documentation of the data dictionary:
- Field descriptions and examples
- Usage guidelines for each field type
- Domain-specific field variations
- Best practices for data entry

### [SCHEMA_VALIDATION_REPORT.md](SCHEMA_VALIDATION_REPORT.md)
Validation report documenting schema conformance:
- Validation results across all domains
- Issues identified and their severity
- Remediation status and recommendations
- Conformance metrics and statistics

## Usage

These schema files are referenced by individual domain schemas located in each resource subdirectory. When adding new resources or modifying existing ones, consult these global definitions to ensure consistency.

## Maintenance

Updates to global schema definitions should be coordinated across all domains to prevent breaking changes. Use the validation report to verify schema conformance after any modifications.
