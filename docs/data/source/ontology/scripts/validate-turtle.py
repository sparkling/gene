#!/usr/bin/env python3
"""
Basic Turtle RDF syntax and conformance validator.
Checks structure without external dependencies.
"""

import re
import sys
from pathlib import Path
from collections import defaultdict

# Validation results
errors = []
warnings = []
stats = defaultdict(int)

def parse_turtle_lite(content: str) -> dict:
    """
    Lightweight Turtle parser for validation.
    Returns dict of subjects with their predicates/objects.
    """
    resources = {}
    current_subject = None
    current_predicates = {}

    # Remove comments
    lines = []
    for line in content.split('\n'):
        # Remove line comments
        if '#' in line:
            # Check if # is inside quotes
            in_string = False
            clean_line = ""
            for i, char in enumerate(line):
                if char == '"' and (i == 0 or line[i-1] != '\\'):
                    in_string = not in_string
                if char == '#' and not in_string:
                    clean_line = line[:i]
                    break
                clean_line = line
            lines.append(clean_line)
        else:
            lines.append(line)

    content = '\n'.join(lines)

    # Extract prefixes
    prefixes = {}
    for match in re.finditer(r'@prefix\s+(\w*:)\s+<([^>]+)>\s*\.', content):
        prefixes[match.group(1)] = match.group(2)

    # Simple resource extraction
    # Match subject blocks ending with .
    resource_pattern = r'(\S+)\s+a\s+(\S+)\s*;([^.]+)\.'
    for match in re.finditer(resource_pattern, content, re.DOTALL):
        subject = match.group(1)
        rdf_type = match.group(2)
        predicates_text = match.group(3)

        resources[subject] = {
            'type': rdf_type,
            'predicates': {}
        }

        # Extract predicates
        for pred_match in re.finditer(r'(\S+)\s+(".*?(?<!\\)"|<[^>]+>|\S+)\s*[;.]?', predicates_text):
            pred = pred_match.group(1)
            obj = pred_match.group(2)
            if pred not in resources[subject]['predicates']:
                resources[subject]['predicates'][pred] = []
            resources[subject]['predicates'][pred].append(obj)

    return {'prefixes': prefixes, 'resources': resources}


def validate_datasource(uri: str, data: dict) -> list:
    """Validate a DataSource instance against SHACL-like constraints."""
    issues = []
    preds = data.get('predicates', {})

    # Required: ds:id
    if 'ds:id' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:id")
    else:
        id_val = preds['ds:id'][0].strip('"')
        if len(id_val) < 2:
            issues.append(f"VIOLATION: {uri} ds:id too short (min 2 chars)")
        if not re.match(r'^[a-z0-9][a-z0-9.-]*[a-z0-9]$', id_val) and len(id_val) > 2:
            issues.append(f"WARNING: {uri} ds:id may not match pattern: {id_val}")

    # Required: ds:title
    if 'ds:title' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:title")

    # Required: ds:belongsToCategory
    if 'ds:belongsToCategory' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:belongsToCategory")

    # Required: ds:belongsToSubcategory
    if 'ds:belongsToSubcategory' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:belongsToSubcategory")

    # Required: ds:status
    if 'ds:status' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:status")
    else:
        status_val = preds['ds:status'][0]
        valid_statuses = ['tax:Draft', 'tax:Active', 'tax:Deprecated']
        if status_val not in valid_statuses:
            issues.append(f"VIOLATION: {uri} invalid status: {status_val}")

    # Recommended: ds:tier
    if 'ds:tier' not in preds:
        issues.append(f"WARNING: {uri} missing recommended ds:tier")
    else:
        tier_val = preds['ds:tier'][0]
        valid_tiers = ['tax:Tier1', 'tax:Tier2', 'tax:Tier3']
        if tier_val not in valid_tiers:
            issues.append(f"VIOLATION: {uri} invalid tier: {tier_val}")

    return issues


def validate_category(uri: str, data: dict) -> list:
    """Validate a Category instance."""
    issues = []
    preds = data.get('predicates', {})

    if 'ds:id' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:id")
    if 'ds:title' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:title")
    if 'ds:categoryClassification' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:categoryClassification")

    return issues


def validate_subcategory(uri: str, data: dict) -> list:
    """Validate a Subcategory instance."""
    issues = []
    preds = data.get('predicates', {})

    if 'ds:id' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:id")
    if 'ds:title' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:title")
    if 'ds:belongsToCategory' not in preds:
        issues.append(f"VIOLATION: {uri} missing required ds:belongsToCategory")

    return issues


def main():
    """Run validation on all ontology files."""
    base_dir = Path(__file__).parent.parent

    files_to_validate = [
        base_dir / "skos" / "datasource-taxonomy.ttl",
        base_dir / "owl" / "datasource-ontology.ttl",
        base_dir / "shacl" / "datasource-shapes.ttl",
        base_dir / "instances" / "datasources.ttl",
    ]

    all_issues = []

    for ttl_file in files_to_validate:
        if not ttl_file.exists():
            print(f"SKIP: {ttl_file.name} not found")
            continue

        print(f"\nValidating: {ttl_file.name}")
        print("-" * 50)

        content = ttl_file.read_text(encoding='utf-8')

        # Basic syntax checks
        open_braces = content.count('[')
        close_braces = content.count(']')
        if open_braces != close_braces:
            all_issues.append(f"SYNTAX: {ttl_file.name} - Mismatched brackets: [ = {open_braces}, ] = {close_braces}")

        # Check for prefix declarations
        prefix_count = len(re.findall(r'@prefix', content))
        print(f"  Prefixes declared: {prefix_count}")

        # Count resources by type
        type_counts = defaultdict(int)
        for match in re.finditer(r'(\S+)\s+a\s+(\S+)\s*[;.]', content):
            rdf_type = match.group(2)
            type_counts[rdf_type] += 1

        print(f"  Resources by type:")
        for rdf_type, count in sorted(type_counts.items()):
            print(f"    {rdf_type}: {count}")

        # Validate instances file
        if "instances" in str(ttl_file):
            parsed = parse_turtle_lite(content)
            resources = parsed.get('resources', {})

            violations = 0
            warnings_count = 0

            for uri, data in resources.items():
                rdf_type = data.get('type', '')

                if rdf_type == 'ds:DataSource':
                    issues = validate_datasource(uri, data)
                elif rdf_type == 'ds:Category':
                    issues = validate_category(uri, data)
                elif rdf_type == 'ds:Subcategory':
                    issues = validate_subcategory(uri, data)
                else:
                    issues = []

                for issue in issues:
                    if issue.startswith("VIOLATION"):
                        violations += 1
                    elif issue.startswith("WARNING"):
                        warnings_count += 1
                    all_issues.append(f"{ttl_file.name}: {issue}")

            print(f"  SHACL-like validation:")
            print(f"    Violations: {violations}")
            print(f"    Warnings: {warnings_count}")

    # Summary
    print("\n" + "=" * 50)
    print("VALIDATION SUMMARY")
    print("=" * 50)

    violation_count = sum(1 for i in all_issues if "VIOLATION" in i or "SYNTAX" in i)
    warning_count = sum(1 for i in all_issues if "WARNING" in i)

    if violation_count == 0 and warning_count == 0:
        print("✓ All files passed validation!")
    else:
        print(f"Violations: {violation_count}")
        print(f"Warnings: {warning_count}")

        if violation_count > 0:
            print("\nViolations:")
            for issue in all_issues:
                if "VIOLATION" in issue or "SYNTAX" in issue:
                    print(f"  • {issue}")

        if warning_count > 0 and warning_count <= 10:
            print("\nWarnings (first 10):")
            warning_shown = 0
            for issue in all_issues:
                if "WARNING" in issue and warning_shown < 10:
                    print(f"  • {issue}")
                    warning_shown += 1

    return 0 if violation_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
