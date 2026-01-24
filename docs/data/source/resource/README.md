# Data Source Resources

Comprehensive catalog of biomedical and life sciences data sources organized by domain.

## Categories

| Category | Description | Sources |
|----------|-------------|---------|
| [01.genetics.genomics](./01.genetics.genomics/) | Variant databases, functional prediction, population genetics, pharmacogenomics | ClinVar, gnomAD, dbSNP, PharmGKB |
| [02.compounds.molecules](./02.compounds.molecules/) | Natural products, pharmaceuticals, drug metabolism, chemical ontology | ChEMBL, DrugBank, PubChem, LOTUS |
| [03.diseases.phenotypes](./03.diseases.phenotypes/) | Disease ontologies, phenotype databases, gene-disease associations | OMIM, HPO, DisGeNET, Open Targets |
| [04.pathways.networks](./04.pathways.networks/) | Metabolic/signaling pathways, protein interactions, regulatory networks | KEGG, Reactome, STRING, GO |
| [05.traditional.medicine](./05.traditional.medicine/) | TCM, Ayurveda, Kampo, Western herbal databases | TCMBANK, IMPPAT, NAPRALERT |
| [06.nutrition.food](./06.nutrition.food/) | Food composition, dietary supplements, bioactive compounds, metabolomics | FooDB, HMDB, USDA FoodData |
| [07.proteins.molecular.biology](./07.proteins.molecular.biology/) | Protein sequences, structures, molecular interactions | UniProt, PDB, AlphaFold DB |
| [08.literature.knowledge](./08.literature.knowledge/) | Scientific literature, knowledge bases, identifier mapping | PubMed, Wikidata, OpenAlex |
| [09.microbiome](./09.microbiome/) | Gut microbiome, body site microbiomes, microbe-host interactions | HMP, GMrepo, HOMD |

## Structure

Each category contains:
- `README.md` - Category overview
- `_category-schema.json` - Unified schema for the category
- `_data-dictionary.md` - Field definitions and mappings
- Subdirectories for each subcategory with individual data sources

Each data source directory contains:
- `README.md` - Source documentation
- `{source}-CLAUDE.md` - LLM-optimized reference (where available)
- Additional technical documentation

## CLAUDE.md Format

LLM-optimized reference files (~50-100 lines) providing:
- Quick metadata (URL, license, update frequency)
- Key identifiers for cross-referencing
- Core schema (essential fields only)
- API endpoints and download URLs
- Practical query examples
- Integration notes

See [templates](../templates/) for the CLAUDE.md template and guide.

## Example CLAUDE.md Files

| File | Data Source | Category |
|------|-------------|----------|
| [chembl-CLAUDE.md](./02.compounds.molecules/2.7.compound.target.interactions/chembl-CLAUDE.md) | ChEMBL | Bioactivity |
| [gnomad-CLAUDE.md](./01.genetics.genomics/1.3.population.genetics/gnomad-CLAUDE.md) | gnomAD | Population Genetics |
| [clinvar-CLAUDE.md](./01.genetics.genomics/1.1.variant.repositories/clinvar-CLAUDE.md) | ClinVar | Clinical Variants |
