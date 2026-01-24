---
id: download-vmh
title: "VMH Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# VMH (Virtual Metabolic Human) Download Instructions

## Quick Start

```bash
# Download human metabolic model (Recon3D)
wget https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D.mat

# Download AGORA2 microbial models
wget https://www.vmh.life/files/reconstructions/AGORA2/AGORA2.tar.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **MATLAB** or **Python** (COBRApy) for model analysis
- **tar** for archive extraction
- 5-50GB storage for all models

## No Registration Required

VMH data is freely available under CC BY 4.0 license.

## Download Methods

### Method 1: Web Interface Downloads

```bash
# Access download page
# Visit https://www.vmh.life/#download

# Human models
wget https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D.mat
wget https://www.vmh.life/files/reconstructions/Recon/3D.01/Recon3D.xml

# AGORA2 complete archive
wget https://www.vmh.life/files/reconstructions/AGORA2/AGORA2.tar.gz
tar -xzf AGORA2.tar.gz
```

### Method 2: Individual Microbe Models

```bash
# Download specific microbe model
# Format: https://www.vmh.life/files/reconstructions/AGORA2/sbml/{Genus_species_strain}.xml

# Example: Bacteroides thetaiotaomicron
wget "https://www.vmh.life/files/reconstructions/AGORA2/sbml/Bacteroides_thetaiotaomicron_VPI_5482.xml"

# Example: Faecalibacterium prausnitzii
wget "https://www.vmh.life/files/reconstructions/AGORA2/sbml/Faecalibacterium_prausnitzii_A2_165.xml"
```

### Method 3: Batch Download AGORA Models

```python
import requests
import os

# Get list of all AGORA models
model_list_url = "https://www.vmh.life/api/microbes"
response = requests.get(model_list_url)
microbes = response.json()

# Download each model
output_dir = "agora_models"
os.makedirs(output_dir, exist_ok=True)

for microbe in microbes[:10]:  # Limit for example
    model_id = microbe['model_id']
    url = f"https://www.vmh.life/files/reconstructions/AGORA2/sbml/{model_id}.xml"

    response = requests.get(url)
    if response.status_code == 200:
        filepath = os.path.join(output_dir, f"{model_id}.xml")
        with open(filepath, 'w') as f:
            f.write(response.text)
        print(f"Downloaded {model_id}")
```

### Method 4: VMH API Access

```bash
# Get reaction information
curl "https://www.vmh.life/api/reactions?page=1&pageSize=100" > reactions.json

# Get metabolite information
curl "https://www.vmh.life/api/metabolites?page=1&pageSize=100" > metabolites.json

# Get microbe list
curl "https://www.vmh.life/api/microbes" > microbes.json

# Get specific microbe model info
curl "https://www.vmh.life/api/microbes/Bacteroides_thetaiotaomicron_VPI_5482" > bt_info.json
```

### Method 5: Using COBRApy

```python
import cobra

# Load model from SBML
model = cobra.io.read_sbml_model("Recon3D.xml")

print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes: {len(model.genes)}")

# Run FBA
solution = model.optimize()
print(f"Growth rate: {solution.objective_value}")

# Save in different format
cobra.io.save_json_model(model, "Recon3D.json")
cobra.io.save_matlab_model(model, "Recon3D_new.mat")
```

### Method 6: Supplementary Tables

```bash
# Download reaction, metabolite, gene tables
wget "https://www.vmh.life/files/tables/reactions.tsv"
wget "https://www.vmh.life/files/tables/metabolites.tsv"
wget "https://www.vmh.life/files/tables/genes.tsv"
wget "https://www.vmh.life/files/tables/microbe_phenotypes.tsv"
```

## File Inventory

### Human Models

| File | Size | Description |
|------|------|-------------|
| Recon3D.mat | ~50 MB | MATLAB format |
| Recon3D.xml | ~100 MB | SBML format |
| Recon3D.json | ~80 MB | JSON format |

### AGORA2 Models

| Archive | Size | Description |
|---------|------|-------------|
| AGORA2.tar.gz | ~5 GB | All 7,206 models |
| Individual models | ~1-5 MB each | Per species |

### Tabular Data

| File | Size | Description |
|------|------|-------------|
| reactions.tsv | ~10 MB | All reactions |
| metabolites.tsv | ~5 MB | All metabolites |
| microbe_phenotypes.tsv | ~2 MB | Phenotype data |

## Post-Download Processing

```python
import cobra
import pandas as pd
import os

# Load and analyze AGORA model
model = cobra.io.read_sbml_model("Bacteroides_thetaiotaomicron_VPI_5482.xml")

# Basic statistics
print(f"Model: {model.id}")
print(f"Reactions: {len(model.reactions)}")
print(f"Metabolites: {len(model.metabolites)}")
print(f"Genes: {len(model.genes)}")

# List exchange reactions (diet/secretion)
exchanges = [r for r in model.reactions if r.id.startswith('EX_')]
print(f"Exchange reactions: {len(exchanges)}")

# Find subsystems
subsystems = set()
for r in model.reactions:
    if r.subsystem:
        subsystems.add(r.subsystem)
print(f"Subsystems: {len(subsystems)}")

# Run FBA
solution = model.optimize()
print(f"Optimal growth: {solution.objective_value:.4f}")

# Get flux distribution
fluxes = solution.fluxes
active_fluxes = fluxes[fluxes.abs() > 1e-6]
print(f"Active reactions: {len(active_fluxes)}")
```

### Create Community Model

```python
import cobra
from cobra import Model, Reaction, Metabolite

def create_community_model(model_files, model_names):
    """Create a simple community model from individual models."""

    community = Model('community')

    for model_file, name in zip(model_files, model_names):
        # Load individual model
        model = cobra.io.read_sbml_model(model_file)

        # Rename reactions and metabolites to avoid conflicts
        for reaction in model.reactions:
            reaction.id = f"{name}_{reaction.id}"
        for metabolite in model.metabolites:
            metabolite.id = f"{name}_{metabolite.id}"

        # Add to community
        community.add_reactions(model.reactions)

    return community

# Example: Create community
models = ["Bacteroides_thetaiotaomicron.xml", "Faecalibacterium_prausnitzii.xml"]
names = ["Bt", "Fp"]
community = create_community_model(models, names)
print(f"Community model: {len(community.reactions)} reactions")
```

### Extract Metabolite Data

```python
import pandas as pd

# Load metabolites table
metabolites = pd.read_csv("metabolites.tsv", sep="\t")

# Filter by criteria
scfa = metabolites[metabolites['name'].str.contains('acetate|butyrate|propionate', case=False)]
print(f"Short-chain fatty acids: {len(scfa)}")

# Get cross-references
kegg_ids = metabolites[['vmh_id', 'kegg_id']].dropna()
print(f"Metabolites with KEGG IDs: {len(kegg_ids)}")
```

## Verification

```bash
# Check SBML validity
python3 << 'EOF'
import cobra
model = cobra.io.read_sbml_model("Recon3D.xml")
print(f"Model loaded successfully: {model.id}")
print(f"Mass balance: {len(cobra.manipulation.validation.check_mass_balance(model))} unbalanced reactions")
EOF

# Verify model can solve
python3 << 'EOF'
import cobra
model = cobra.io.read_sbml_model("Recon3D.xml")
solution = model.optimize()
print(f"Optimization status: {solution.status}")
print(f"Objective value: {solution.objective_value}")
EOF
```

## Update Schedule

| Release | Frequency |
|---------|-----------|
| AGORA updates | Major versions (AGORA2, etc.) |
| Recon updates | Major versions (Recon3D, etc.) |
| Bug fixes | As needed |

## Common Issues

- **Large files**: Use chunked downloads for AGORA archive
- **SBML parsing**: Ensure COBRApy is up to date
- **Solver requirements**: Install glpk, gurobi, or cplex
- **Memory usage**: Large community models need significant RAM
- **Model compatibility**: Check SBML level/version

## Software Requirements

```bash
# Install COBRApy
pip install cobra

# Install solver (glpk is open source)
conda install -c conda-forge glpk

# Or use Gurobi (academic license available)
# Or use CPLEX (IBM academic initiative)
```

## Integration Examples

```python
# Link to gut microbiome data
import cobra
import pandas as pd

# Load microbiome abundances from GMrepo or HMP
abundances = pd.read_csv("sample_abundances.tsv", sep="\t", index_col=0)

# Map species to AGORA models
agora_mapping = pd.read_csv("agora_taxonomy.tsv", sep="\t")

# Create personalized community model
def build_personalized_model(abundances, top_n=20):
    """Build community model from top N abundant species."""
    top_species = abundances.nlargest(top_n)

    models = []
    for species, abundance in top_species.items():
        # Find matching AGORA model
        agora_match = agora_mapping[
            agora_mapping['species'].str.contains(species, case=False)
        ]
        if len(agora_match) > 0:
            model_id = agora_match.iloc[0]['model_id']
            models.append((model_id, abundance))

    return models

personalized = build_personalized_model(abundances.iloc[:, 0])
print(f"Selected {len(personalized)} models for personalized analysis")
```

## Related Resources

- [KEGG](../../../../04.pathways.networks/4.1.metabolic.pathways/kegg/_index.md) - Pathway data
- [Reactome](../../../../04.pathways.networks/4.1.metabolic.pathways/reactome/_index.md) - Reaction database
- [GMrepo](../../9.1.gut.microbiome/gmrepo/_index.md) - Microbiome abundances
- [gutMDisorder](../gutmdisorder/_index.md) - Disease associations
