#!/usr/bin/env python3
"""
Generate RDF instance data from data source README frontmatter.
Conforms to datasource-ontology.ttl and datasource-shapes.ttl.
"""

import os
import re
import yaml
from pathlib import Path
from datetime import date

# Paths
RESOURCE_DIR = Path(__file__).parent.parent.parent / "resource"
OUTPUT_FILE = Path(__file__).parent.parent / "instances" / "datasources.ttl"

# Namespace mappings for categories
CATEGORY_MAP = {
    "genetics.genomics": "GeneticsGenomics",
    "compounds.molecules": "CompoundsMolecules",
    "diseases.phenotypes": "DiseasesPhenotypes",
    "pathways.networks": "PathwaysNetworks",
    "traditional.medicine": "TraditionalMedicine",
    "nutrition.food": "NutritionFood",
    "proteins.molecular.biology": "ProteinsMolecularBiology",
    "literature.knowledge": "LiteratureKnowledge",
    "microbiome": "Microbiome",
}

# Subcategory mappings
SUBCATEGORY_MAP = {
    "variant.repositories": "VariantRepositories",
    "functional.prediction": "FunctionalPrediction",
    "population.genetics": "PopulationGenetics",
    "pharmacogenomics": "Pharmacogenomics",
    "expression.regulation": "ExpressionRegulation",
    "cancer.genomics": "CancerGenomics",
    "natural.products": "NaturalProducts",
    "pharmaceuticals": "Pharmaceuticals",
    "traditional.medicine.compounds": "TraditionalMedicineCompounds",
    "food.compounds.nutrients": "FoodCompoundsNutrients",
    "drug.metabolism.pharmacokinetics": "DrugMetabolism",
    "chemical.ontology.classification": "ChemicalOntology",
    "compound.target.interactions": "CompoundTargetInteractions",
    "disease.ontologies": "DiseaseOntologies",
    "phenotype.databases": "PhenotypeDatabases",
    "disease.gene.associations": "DiseaseGeneAssociations",
    "cancer.oncology": "CancerOncology",
    "rare.orphan.diseases": "RareOrphanDiseases",
    "autoimmune.inflammatory": "AutoimmuneInflammatory",
    "mental.health.neurological": "MentalHealthNeurological",
    "metabolic.pathways": "MetabolicPathways",
    "signaling.pathways": "SignalingPathways",
    "protein.protein.interactions": "ProteinProteinInteractions",
    "drug.target.interactions": "DrugTargetInteractions",
    "gene.function.ontology": "GeneFunctionOntology",
    "regulatory.networks": "RegulatoryNetworks",
    "traditional.chinese.medicine": "TraditionalChineseMedicine",
    "south.east.asian.systems": "SouthEastAsianSystems",
    "western.global.herbal": "WesternGlobalHerbal",
    "multi.system.integration": "MultiSystemIntegration",
    "food.composition": "FoodComposition",
    "dietary.supplements": "DietarySupplements",
    "bioactive.food.compounds": "BioactiveFoodCompounds",
    "metabolomics": "Metabolomics",
    "protein.sequences.annotations": "ProteinSequencesAnnotations",
    "protein.structures": "ProteinStructures",
    "molecular.interactions": "MolecularInteractions",
    "scientific.literature": "ScientificLiterature",
    "knowledge.bases": "KnowledgeBases",
    "identifier.mapping": "IdentifierMapping",
    "regulatory.legal": "RegulatoryLegal",
    "gut.microbiome": "GutMicrobiome",
    "body.site.microbiomes": "BodySiteMicrobiomes",
    "microbe.host.interactions": "MicrobeHostInteractions",
}

# Tier mapping
TIER_MAP = {
    1: "Tier1",
    "1": "Tier1",
    2: "Tier2",
    "2": "Tier2",
    3: "Tier3",
    "3": "Tier3",
}

# Status mapping
STATUS_MAP = {
    "active": "Active",
    "draft": "Draft",
    "deprecated": "Deprecated",
}


def extract_frontmatter(readme_path: Path) -> dict:
    """Extract YAML frontmatter from README.md."""
    content = readme_path.read_text(encoding="utf-8")
    match = re.match(r"^---\n(.*?)\n---", content, re.DOTALL)
    if match:
        try:
            return yaml.safe_load(match.group(1))
        except yaml.YAMLError:
            return {}
    return {}


def escape_turtle_string(s: str) -> str:
    """Escape special characters for Turtle string literals."""
    if not s:
        return ""
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")


def generate_source_turtle(data: dict, source_path: Path) -> str:
    """Generate Turtle RDF for a data source."""
    source_id = data.get("id", "")
    if not source_id:
        return ""

    # Clean the ID for use as URI local name
    uri_id = re.sub(r"[^a-zA-Z0-9_-]", "-", source_id)

    title = escape_turtle_string(data.get("title", source_id))
    category = data.get("category", "")
    subcategory = data.get("subcategory", "")
    tier = data.get("tier", 2)
    status = data.get("status", "draft")
    tags = data.get("tags", [])
    last_updated = data.get("last_updated", str(date.today()))

    # Map to taxonomy concepts
    category_concept = CATEGORY_MAP.get(category, "")
    subcategory_concept = SUBCATEGORY_MAP.get(subcategory, "")
    tier_concept = TIER_MAP.get(tier, "Tier2")
    status_concept = STATUS_MAP.get(status, "Draft")

    if not category_concept or not subcategory_concept:
        print(f"Warning: Unknown category/subcategory for {source_id}: {category}/{subcategory}")
        return ""

    lines = [
        f"data:{uri_id} a ds:DataSource ;",
        f'    ds:id "{source_id}" ;',
        f'    ds:title "{title}" ;',
        f"    ds:belongsToCategory cat:{category_concept} ;",
        f"    ds:belongsToSubcategory subcat:{subcategory_concept} ;",
        f"    ds:tier tax:{tier_concept} ;",
        f"    ds:status tax:{status_concept} ;",
    ]

    # Add tags
    for tag in tags:
        clean_tag = escape_turtle_string(tag)
        lines.append(f'    ds:tag "{clean_tag}" ;')

    # Add last updated
    if last_updated:
        lines.append(f'    ds:lastUpdated "{last_updated}"^^xsd:date ;')

    # Add source path as description
    rel_path = source_path.relative_to(RESOURCE_DIR)
    lines.append(f'    rdfs:comment "Source documentation: {rel_path}" ;')

    # Close the resource
    lines[-1] = lines[-1].rstrip(" ;") + " ."

    return "\n".join(lines)


def generate_category_instances() -> str:
    """Generate category and subcategory instances."""
    lines = []

    # Categories
    for cat_key, cat_concept in CATEGORY_MAP.items():
        lines.append(f"cat:{cat_concept} a ds:Category ;")
        lines.append(f'    ds:id "{cat_key}" ;')
        lines.append(f'    ds:title "{cat_concept}" ;')
        lines.append(f"    ds:categoryClassification tax:{cat_concept} ;")
        lines.append(f"    ds:status tax:Active .")
        lines.append("")

    # Subcategories
    for subcat_key, subcat_concept in SUBCATEGORY_MAP.items():
        # Find parent category
        parent_cat = None
        for cat_key, cat_concept in CATEGORY_MAP.items():
            # Check if this subcategory belongs to this category
            # This is a heuristic based on folder structure
            if subcat_key in ["variant.repositories", "functional.prediction", "population.genetics",
                             "pharmacogenomics", "expression.regulation", "cancer.genomics"]:
                parent_cat = "GeneticsGenomics"
                break
            elif subcat_key in ["natural.products", "pharmaceuticals", "traditional.medicine.compounds",
                               "food.compounds.nutrients", "drug.metabolism.pharmacokinetics",
                               "chemical.ontology.classification", "compound.target.interactions"]:
                parent_cat = "CompoundsMolecules"
                break
            elif subcat_key in ["disease.ontologies", "phenotype.databases", "disease.gene.associations",
                               "cancer.oncology", "rare.orphan.diseases", "autoimmune.inflammatory",
                               "mental.health.neurological"]:
                parent_cat = "DiseasesPhenotypes"
                break
            elif subcat_key in ["metabolic.pathways", "signaling.pathways", "protein.protein.interactions",
                               "drug.target.interactions", "gene.function.ontology", "regulatory.networks"]:
                parent_cat = "PathwaysNetworks"
                break
            elif subcat_key in ["traditional.chinese.medicine", "south.east.asian.systems",
                               "western.global.herbal", "multi.system.integration"]:
                parent_cat = "TraditionalMedicine"
                break
            elif subcat_key in ["food.composition", "dietary.supplements", "bioactive.food.compounds",
                               "metabolomics"]:
                parent_cat = "NutritionFood"
                break
            elif subcat_key in ["protein.sequences.annotations", "protein.structures", "molecular.interactions"]:
                parent_cat = "ProteinsMolecularBiology"
                break
            elif subcat_key in ["scientific.literature", "knowledge.bases", "identifier.mapping",
                               "regulatory.legal"]:
                parent_cat = "LiteratureKnowledge"
                break
            elif subcat_key in ["gut.microbiome", "body.site.microbiomes", "microbe.host.interactions"]:
                parent_cat = "Microbiome"
                break

        if not parent_cat:
            continue

        lines.append(f"subcat:{subcat_concept} a ds:Subcategory ;")
        lines.append(f'    ds:id "{subcat_key}" ;')
        lines.append(f'    ds:title "{subcat_concept}" ;')
        lines.append(f"    ds:belongsToCategory cat:{parent_cat} ;")
        lines.append(f"    ds:subcategoryClassification tax:{subcat_concept} ;")
        lines.append(f"    ds:status tax:Active .")
        lines.append("")

    return "\n".join(lines)


def main():
    """Generate instance data from all README files."""

    # Collect all data sources
    sources = []
    for readme_path in RESOURCE_DIR.rglob("README.md"):
        # Only process source-level READMEs (depth 3)
        rel_path = readme_path.relative_to(RESOURCE_DIR)
        parts = rel_path.parts
        if len(parts) == 4:  # category/subcategory/source/README.md
            data = extract_frontmatter(readme_path)
            if data.get("type") == "source" and data.get("id"):
                sources.append((data, readme_path.parent))

    print(f"Found {len(sources)} data sources")

    # Generate Turtle output
    output_lines = [
        "# =============================================================================",
        "# Gene Data Source Instance Data",
        "# =============================================================================",
        f"# Generated: {date.today()}",
        f"# Sources: {len(sources)}",
        "# Conforms to: datasource-ontology.ttl, datasource-shapes.ttl",
        "# =============================================================================",
        "",
        "@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .",
        "@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .",
        "@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .",
        "@prefix ds: <https://gene.ai/ontology/datasource#> .",
        "@prefix tax: <https://gene.ai/taxonomy/> .",
        "@prefix data: <https://gene.ai/data/source/> .",
        "@prefix cat: <https://gene.ai/data/category/> .",
        "@prefix subcat: <https://gene.ai/data/subcategory/> .",
        "",
        "# =============================================================================",
        "# CATEGORIES",
        "# =============================================================================",
        "",
    ]

    # Add category and subcategory instances
    output_lines.append(generate_category_instances())

    output_lines.extend([
        "",
        "# =============================================================================",
        "# DATA SOURCES",
        "# =============================================================================",
        "",
    ])

    # Add data source instances
    for data, source_path in sorted(sources, key=lambda x: x[0].get("id", "")):
        turtle = generate_source_turtle(data, source_path)
        if turtle:
            output_lines.append(turtle)
            output_lines.append("")

    # Write output
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text("\n".join(output_lines), encoding="utf-8")
    print(f"Generated: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
