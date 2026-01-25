---
id: stitch
title: "STITCH - Chemical-Protein Interaction Networks"
type: source
parent: ../README.md
tier: 2
status: active
category: pathways.networks
subcategory: drug.target.interactions
tags:
  - drug-target
  - chemical-protein
  - interactions
  - network
  - open-access
---

# STITCH - Chemical-Protein Interaction Networks

**Category:** [Pathways & Networks](../../README.md) > [Drug-Target Interactions](../README.md)

## Overview

STITCH (Search Tool for Interacting Chemicals) is a database of known and predicted interactions between chemicals and proteins. It is the chemical-protein counterpart to STRING, using similar evidence integration methods to compute interaction confidence scores. STITCH combines experimental data, pathway databases, text mining, and chemical structure similarity to predict which proteins a given chemical may interact with.

The database integrates information from metabolic pathways, crystal structures, binding experiments, drug-target databases, and literature mining. STITCH covers over 500,000 chemicals and 9.6 million proteins across 2,031 organisms, with 1.6 billion chemical-protein interactions.

STITCH is maintained by the same team that develops STRING and uses compatible identifiers and scoring systems, enabling seamless integration of chemical-protein and protein-protein interaction networks.

## Key Statistics

| Metric | Value |
|--------|-------|
| Chemicals | 500,000+ |
| Proteins | 9,643,763 |
| Organisms | 2,031 |
| Chemical-Protein Interactions | 1.6+ billion |
| Drug-like Chemicals | 350,000+ |
| Metabolites | 150,000+ |

## Primary Use Cases

1. **Drug target discovery** - Identify potential protein targets for compounds
2. **Off-target prediction** - Find potential side-effect mechanisms
3. **Polypharmacology analysis** - Study drugs with multiple targets
4. **Drug repurposing** - Find new indications via target networks
5. **Metabolite-protein networks** - Study metabolic regulation

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| STITCH Chemical | `CID{m/s}{#########}` | CIDm00002244 |
| STRING Protein | `{taxid}.{ENSP}` | 9606.ENSP00000269305 |
| PubChem CID | Numeric | 2244 |
| DrugBank | `DB{#####}` | DB00945 |
| ChEBI | `CHEBI:{#####}` | CHEBI:15365 |

### Chemical ID Prefixes

| Prefix | Meaning |
|--------|---------|
| CIDm | Merged/stereo-merged compound |
| CIDs | Stereospecific compound |

## Score Channels

| Score | Channel | Description |
|-------|---------|-------------|
| experimental | Experimental | Binding assays, crystal structures |
| database | Database | DrugBank, KEGG, pathway databases |
| textmining | Text Mining | Literature co-occurrence |
| prediction | Prediction | Binding site similarity, docking |
| combined | Combined | Probabilistic integration |

### Score Confidence Levels

| Combined Score | Confidence |
|---------------|------------|
| 0.900 - 1.000 | Highest |
| 0.700 - 0.899 | High |
| 0.400 - 0.699 | Medium |
| 0.150 - 0.399 | Low |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | http://stitch.embl.de | Interactive networks |
| REST API | http://stitch.embl.de/api/ | JSON/TSV/XML |
| Downloads | http://stitch.embl.de/download/ | Bulk data |

### API Examples

```bash
# Get interactions for a chemical
curl "http://stitch.embl.de/api/json/interactions?identifier=CIDm00002244&species=9606&limit=10"

# Get interactions for aspirin by name
curl "http://stitch.embl.de/api/json/interactions?identifier=aspirin&species=9606"

# Map chemical name to STITCH ID
curl "http://stitch.embl.de/api/json/resolve?identifier=aspirin"

# Get network with proteins
curl "http://stitch.embl.de/api/json/network?identifiers=CIDm00002244%0dTP53&species=9606"

# Functional enrichment for targets
curl "http://stitch.embl.de/api/json/enrichment?identifiers=CIDm00002244&species=9606"

# Get network image
curl "http://stitch.embl.de/api/image/network?identifier=aspirin&species=9606" > network.png
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| JSON | Structured | API responses |
| TSV | Tab-delimited | Bulk processing |
| XML | Structured | Legacy integration |
| PNG/SVG | Images | Visualization |

## API Response Fields

### Interaction Response

| Field | Type | Description |
|-------|------|-------------|
| chemicalId | String | STITCH chemical ID |
| stringId | String | STRING protein ID |
| preferredName_chemical | String | Chemical name |
| preferredName_protein | String | Gene symbol |
| ncbiTaxonId | Integer | Taxonomy ID |
| score | Float | Combined score |
| experimentalScore | Float | Experimental evidence |
| databaseScore | Float | Database evidence |
| textminingScore | Float | Text mining evidence |

## Data Sources

| Source | Type | Content |
|--------|------|---------|
| ChEMBL | Bioactivity | Binding assays |
| DrugBank | Drug-target | Approved drugs |
| KEGG | Pathway | Metabolic targets |
| PDB | Structure | Crystal structures |
| CTD | Toxicology | Chemical-gene associations |
| Matador | Drug-target | Drug targets |
| PubChem | Bioassay | HTS data |
| PDSP | Binding | Ki database |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |

## Cross-References

| Database | Relationship |
|----------|--------------|
| PubChem | Chemical identifiers |
| ChEBI | Chemical ontology |
| DrugBank | Drug information |
| KEGG | Compound IDs |
| STRING | Protein IDs |
| UniProt | Protein accessions |
| Ensembl | Gene IDs |

## Integration with STRING

STITCH and STRING use compatible systems:

- Same protein identifiers (taxid.ENSP)
- Same confidence scoring framework
- Combined chemical-protein-protein networks
- Shared API structure

```bash
# Combined network query (STITCH + STRING)
curl "http://stitch.embl.de/api/json/network?identifiers=aspirin%0dTP53%0dMDM2&species=9606"
```

## Bulk Download Files

| File | Description |
|------|-------------|
| chemical_chemical.links.v5.0.tsv.gz | Chemical-chemical interactions |
| protein_chemical.links.v5.0.tsv.gz | Protein-chemical interactions |
| chemicals.v5.0.tsv.gz | Chemical annotations |
| chemical.sources.v5.0.tsv.gz | Evidence sources |

## Limitations

- HTTP (not HTTPS) may cause security issues in some environments
- Predicted interactions may include false positives
- Chemical ID mapping can be complex (merged vs stereospecific)
- Text-mining scores may not reflect direct binding evidence

## Actions (Interaction Types)

| Action | Description |
|--------|-------------|
| activation | Activates protein |
| inhibition | Inhibits protein |
| binding | Binds to protein |
| catalysis | Substrate/product |
| reaction | Metabolic reaction |
| expression | Affects expression |

## See Also

- [Schema Documentation](./schema.md) - Technical schema and data formats
- [Download Instructions](./download.md) - Bulk data acquisition guide
- [STRING](../../4.3.protein.protein.interactions/string/README.md) - Protein-protein interactions
- [ChEMBL](../../../03.bioactivity.drug.data/3.1.compound.bioactivity/chembl/README.md) - Bioactivity data
- [DrugBank](../../../03.bioactivity.drug.data/3.2.drug.databases/drugbank/README.md) - Drug information
