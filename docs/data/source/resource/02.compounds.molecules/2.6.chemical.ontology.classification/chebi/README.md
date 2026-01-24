---
id: chebi
title: "ChEBI - Chemical Entities of Biological Interest"
type: data-source
category: compounds
subcategory: chemical-ontology
parent: ../_index.md
tier: 1
last_updated: 2026-01-23
status: active
tags: [ontology, chemical-classification, biological-roles, molecular-entities, obo]
---

# ChEBI - Chemical Entities of Biological Interest

**Category:** [Compounds & Molecules](../../_index.md) > [Chemical Ontology & Classification](../_index.md)

## Overview

ChEBI (Chemical Entities of Biological Interest) is a freely available dictionary and ontology of molecular entities focused on small chemical compounds relevant to biological systems. Developed and maintained by the European Bioinformatics Institute (EMBL-EBI), ChEBI provides standardized nomenclature and classification for chemical entities used in life science databases.

Unlike simple chemical databases, ChEBI classifies compounds according to their biological roles (e.g., enzyme inhibitor, hormone), chemical structure (e.g., steroid, alkaloid), and relationships to other entities (e.g., conjugate acid/base, enantiomer). This ontological structure enables sophisticated queries about chemical relationships and biological functions.

ChEBI is a key component of the Open Biomedical Ontologies (OBO) Foundry and integrates with other biological ontologies including the Gene Ontology, enabling cross-ontology queries linking chemistry to biology.

## Key Statistics

| Metric | Value |
|--------|-------|
| Chemical Entities | 150,000+ |
| 3-Star Manually Curated | 61,000+ |
| Relationships | 500,000+ |
| Ontology Classes | 6,000+ |
| Synonyms | 600,000+ |
| Cross-References | 100,000+ |
| Updates | Monthly |

## Primary Use Cases

1. **Chemical Annotation** - Standardize chemical entity references in databases
2. **Ontology Queries** - Find chemicals by role, structure class, or relationship
3. **Data Integration** - Link chemical data across biological databases
4. **Text Mining** - Named entity recognition for chemicals in literature
5. **Pathway Annotation** - Chemical entity classification for metabolic pathways

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| ChEBI ID | CHEBI: + digits | CHEBI:15377 (water) |
| InChI Key | 27 characters | Standard format |
| Name | Text | Standardized nomenclature |
| SMILES | Variable | Structure representation |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.ebi.ac.uk/chebi/ | Search and browse |
| Ontology Files | https://www.ebi.ac.uk/chebi/downloadsForward.do | OWL, OBO formats |
| REST API | https://www.ebi.ac.uk/webservices/chebi | Programmatic access |
| SPARQL Endpoint | https://www.ebi.ac.uk/rdf/services/sparql | Semantic queries |

## Ontology Relationships

| Relationship | Description |
|--------------|-------------|
| is_a | Subsumption (subclass) |
| has_role | Biological/chemical role |
| has_part | Structural component |
| is_conjugate_acid_of | Acid-base pair |
| is_enantiomer_of | Stereochemical relationship |
| is_tautomer_of | Tautomeric forms |

## Data Quality Stars

| Stars | Description |
|-------|-------------|
| 3 | Fully manually curated |
| 2 | Automatically derived, expert reviewed |
| 1 | Automatically derived |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution Required | Yes |

## Related Resources

- [PubChem](../pubchem/_index.md) - Chemical repository
- [ClassyFire](../classyfire/_index.md) - Chemical taxonomy
- [ChEMBL](../../2.2.pharmaceuticals/chembl/_index.md) - Bioactivity

## See Also

- [Schema Documentation](./schema.md)

## References

1. Hastings J, et al. (2016) "ChEBI in 2016: Improved services and an expanding collection of metabolites." Nucleic Acids Res. 44(D1):D1214-D1219.
