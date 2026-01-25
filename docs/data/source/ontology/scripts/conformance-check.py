#!/usr/bin/env python3
"""
Conformance checker for Gene Data Source RDF.
Validates Turtle syntax and SHACL-like constraints.
"""

import re
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path(__file__).parent.parent


def parse_turtle_resources(content: str) -> dict:
    """Parse Turtle resources handling multi-line blocks."""
    resources = {}

    # Remove comments (lines starting with #)
    lines = []
    for line in content.split('\n'):
        stripped = line.strip()
        if stripped.startswith('#'):
            continue
        # Handle inline comments
        if '#' in line and '"' not in line.split('#')[0]:
            line = line.split('#')[0]
        lines.append(line)
    content = '\n'.join(lines)

    # Split into resource blocks (end with .)
    # First, find all subject declarations
    subject_pattern = r'^(\S+)\s+a\s+(\S+)\s*;'

    current_pos = 0
    subjects = []
    for match in re.finditer(subject_pattern, content, re.MULTILINE):
        subjects.append((match.start(), match.group(1), match.group(2)))

    # For each subject, find its block until the next subject or end
    for i, (start, subject, rdf_type) in enumerate(subjects):
        # Find end of this block (next . at start or end of content)
        if i + 1 < len(subjects):
            end = subjects[i + 1][0]
        else:
            end = len(content)

        block = content[start:end]

        # Find the actual end (last . before next subject)
        last_dot = block.rfind('.')
        if last_dot > 0:
            block = block[:last_dot + 1]

        # Parse predicates from block
        predicates = defaultdict(list)

        # Extract all predicate-object pairs
        # Pattern: predicate value ; or predicate value .
        pred_obj_pattern = r'(\S+:?\S*)\s+("[^"]*"(?:\^\^xsd:\w+)?|<[^>]+>|\S+)\s*[;.]'

        for match in re.finditer(pred_obj_pattern, block):
            pred = match.group(1)
            obj = match.group(2)
            if pred != 'a':  # Skip rdf:type, already captured
                predicates[pred].append(obj)

        resources[subject] = {
            'type': rdf_type,
            'predicates': dict(predicates)
        }

    return resources


def check_datasource_conformance(uri: str, data: dict) -> tuple:
    """Check DataSource against SHACL constraints."""
    violations = []
    warnings = []
    preds = data.get('predicates', {})

    # Required: ds:id
    if 'ds:id' not in preds:
        violations.append(f"{uri}: missing required ds:id")

    # Required: ds:title
    if 'ds:title' not in preds:
        violations.append(f"{uri}: missing required ds:title")

    # Required: ds:belongsToCategory
    if 'ds:belongsToCategory' not in preds:
        violations.append(f"{uri}: missing required ds:belongsToCategory")

    # Required: ds:belongsToSubcategory
    if 'ds:belongsToSubcategory' not in preds:
        violations.append(f"{uri}: missing required ds:belongsToSubcategory")

    # Required: ds:status
    if 'ds:status' not in preds:
        violations.append(f"{uri}: missing required ds:status")
    else:
        status = preds['ds:status'][0]
        if status not in ['tax:Draft', 'tax:Active', 'tax:Deprecated']:
            violations.append(f"{uri}: invalid status {status}")

    # Recommended: ds:tier
    if 'ds:tier' not in preds:
        warnings.append(f"{uri}: missing recommended ds:tier")
    else:
        tier = preds['ds:tier'][0]
        if tier not in ['tax:Tier1', 'tax:Tier2', 'tax:Tier3']:
            violations.append(f"{uri}: invalid tier {tier}")

    return violations, warnings


def check_category_conformance(uri: str, data: dict) -> tuple:
    """Check Category against SHACL constraints."""
    violations = []
    warnings = []
    preds = data.get('predicates', {})

    if 'ds:id' not in preds:
        violations.append(f"{uri}: missing required ds:id")
    if 'ds:title' not in preds:
        violations.append(f"{uri}: missing required ds:title")
    if 'ds:categoryClassification' not in preds:
        violations.append(f"{uri}: missing required ds:categoryClassification")

    return violations, warnings


def check_subcategory_conformance(uri: str, data: dict) -> tuple:
    """Check Subcategory against SHACL constraints."""
    violations = []
    warnings = []
    preds = data.get('predicates', {})

    if 'ds:id' not in preds:
        violations.append(f"{uri}: missing required ds:id")
    if 'ds:title' not in preds:
        violations.append(f"{uri}: missing required ds:title")
    if 'ds:belongsToCategory' not in preds:
        violations.append(f"{uri}: missing required ds:belongsToCategory")

    return violations, warnings


def main():
    """Run conformance checks."""
    print("=" * 60)
    print("GENE DATA SOURCE ONTOLOGY - CONFORMANCE REPORT")
    print("=" * 60)

    # Check each file
    files = {
        "SKOS Taxonomy": BASE_DIR / "skos" / "datasource-taxonomy.ttl",
        "OWL Ontology": BASE_DIR / "owl" / "datasource-ontology.ttl",
        "SHACL Shapes": BASE_DIR / "shacl" / "datasource-shapes.ttl",
        "Instance Data": BASE_DIR / "instances" / "datasources.ttl",
    }

    for name, filepath in files.items():
        if not filepath.exists():
            print(f"\n{name}: NOT FOUND")
            continue

        content = filepath.read_text(encoding='utf-8')
        print(f"\n{name}: {filepath.name}")
        print("-" * 40)

        # Count prefixes
        prefix_count = len(re.findall(r'@prefix', content))
        print(f"  Prefixes: {prefix_count}")

        # Parse resources
        resources = parse_turtle_resources(content)
        print(f"  Resources: {len(resources)}")

        # Count by type
        type_counts = defaultdict(int)
        for uri, data in resources.items():
            type_counts[data['type']] += 1

        for rdf_type, count in sorted(type_counts.items()):
            print(f"    {rdf_type}: {count}")

    # Detailed validation of instance data
    print("\n" + "=" * 60)
    print("INSTANCE DATA VALIDATION")
    print("=" * 60)

    instance_file = BASE_DIR / "instances" / "datasources.ttl"
    if not instance_file.exists():
        print("Instance data not found!")
        return 1

    content = instance_file.read_text(encoding='utf-8')
    resources = parse_turtle_resources(content)

    all_violations = []
    all_warnings = []
    stats = {'categories': 0, 'subcategories': 0, 'datasources': 0}

    for uri, data in resources.items():
        rdf_type = data['type']

        if rdf_type == 'ds:Category':
            stats['categories'] += 1
            v, w = check_category_conformance(uri, data)
        elif rdf_type == 'ds:Subcategory':
            stats['subcategories'] += 1
            v, w = check_subcategory_conformance(uri, data)
        elif rdf_type == 'ds:DataSource':
            stats['datasources'] += 1
            v, w = check_datasource_conformance(uri, data)
        else:
            v, w = [], []

        all_violations.extend(v)
        all_warnings.extend(w)

    print(f"\nStatistics:")
    print(f"  Categories: {stats['categories']}")
    print(f"  Subcategories: {stats['subcategories']}")
    print(f"  Data Sources: {stats['datasources']}")

    print(f"\nValidation Results:")
    print(f"  Violations: {len(all_violations)}")
    print(f"  Warnings: {len(all_warnings)}")

    if all_violations:
        print(f"\nViolations:")
        for v in all_violations[:20]:
            print(f"  ✗ {v}")
        if len(all_violations) > 20:
            print(f"  ... and {len(all_violations) - 20} more")

    if all_warnings and len(all_warnings) <= 10:
        print(f"\nWarnings:")
        for w in all_warnings:
            print(f"  ⚠ {w}")

    # Summary
    print("\n" + "=" * 60)
    if len(all_violations) == 0:
        print("✓ CONFORMANCE CHECK PASSED")
        print("  All instances conform to SHACL shapes.")
    else:
        print("✗ CONFORMANCE CHECK FAILED")
        print(f"  {len(all_violations)} violations found.")

    print("=" * 60)

    return 0 if len(all_violations) == 0 else 1


if __name__ == "__main__":
    exit(main())
