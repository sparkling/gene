# Bulk Download Methods for Pathway and Gene/Protein Target Databases

This document provides comprehensive information on bulk download methods for major biological pathway and gene/protein target databases. Last updated: January 2026.

---

## Table of Contents

1. [Reactome](#1-reactome)
2. [KEGG](#2-kegg)
3. [WikiPathways](#3-wikipathways)
4. [UniProt](#4-uniprot)
5. [STRING](#5-string)
6. [Gene Ontology](#6-gene-ontology)
7. [NCBI Gene](#7-ncbi-gene)
8. [Ensembl](#8-ensembl)
9. [Processing Recommendations Summary](#9-processing-recommendations-summary)

---

## 1. Reactome

**Website:** https://reactome.org/
**Download Page:** https://reactome.org/download-data
**License:** Creative Commons Attribution 4.0

### Core Data Files

| File | URL | Description |
|------|-----|-------------|
| ReactomePathways.txt | https://reactome.org/download/current/ReactomePathways.txt | Complete list of pathways describing the Reactome hierarchy |
| ReactomePathwaysRelation.txt | https://reactome.org/download/current/ReactomePathwaysRelation.txt | Pathways hierarchy relationship documentation |
| NCBI2Reactome.txt | https://reactome.org/download/current/NCBI2Reactome.txt | Maps NCBI Gene IDs to lowest-level pathway diagrams |
| UniProt2Reactome.txt | https://reactome.org/download/current/UniProt2Reactome.txt | Links UniProt accessions to pathway information |
| Ensembl2Reactome.txt | https://reactome.org/download/current/Ensembl2Reactome.txt | Maps Ensembl gene IDs to pathways |
| ChEBI2Reactome.txt | https://reactome.org/download/current/ChEBI2Reactome.txt | Maps ChEBI compound IDs to pathways |

### Graph Database & Advanced Formats

| File | URL | Description | Size |
|------|-----|-------------|------|
| reactome.graphdb.tgz | https://reactome.org/download/current/reactome.graphdb.tgz | Neo4j graph database | ~2-3 GB |
| biopax.zip | https://reactome.org/download/current/biopax.zip | BioPAX Level 3 format | ~500 MB |
| homo_sapiens.3.1.sbml.tgz | https://reactome.org/download/current/homo_sapiens.3.1.sbml.tgz | Human SBML (Level 3, v1) | ~100 MB |

### Update Frequency

- **Quarterly releases** (approximately every 3 months)
- Zenodo repository for archived releases starting from version 89 (June 2024)

### Processing Recommendations

```bash
# Download core mapping files
wget https://reactome.org/download/current/UniProt2Reactome.txt
wget https://reactome.org/download/current/NCBI2Reactome.txt
wget https://reactome.org/download/current/ReactomePathways.txt
wget https://reactome.org/download/current/ReactomePathwaysRelation.txt

# For Neo4j database
wget https://reactome.org/download/current/reactome.graphdb.tgz
tar -xzf reactome.graphdb.tgz
# Import into Neo4j 4.x or 5.x
```

**Best Practices:**
- Use UniProt2Reactome.txt for protein-centric pathway mapping
- Use NCBI2Reactome.txt for gene-centric analysis
- Neo4j graph database recommended for complex queries and relationship traversal
- BioPAX format useful for integration with other pathway tools (e.g., Cytoscape)

---

## 2. KEGG

**Website:** https://www.kegg.jp/
**Download Info:** https://www.kegg.jp/kegg/download/
**License:** Academic subscription required for full FTP access

### Access Tiers

| Tier | Access | Data Available |
|------|--------|----------------|
| Free | REST API | Limited queries (max 10 entries per request) |
| Free | GenomeNet FTP | Subset of KEGG MEDICUS |
| Academic | Subscription FTP | Full KEGG data (weekly updates) |
| Commercial | License Agreement | Full access with commercial use rights |

### KEGG REST API (Free)

**Base URL:** https://rest.kegg.jp/

| Operation | Endpoint | Description | Limitations |
|-----------|----------|-------------|-------------|
| info | /info/{database} | Database release information | - |
| list | /list/{database} | List entry identifiers | Max 10 entries |
| find | /find/{database}/{query} | Search for entries | - |
| get | /get/{entry-ids} | Retrieve entries | Max 10 entries per request |
| conv | /conv/{target}/{source} | Convert identifiers | Specific database pairs only |
| link | /link/{target}/{source} | Cross-references | - |
| ddi | /ddi/{drug-ids} | Drug-drug interactions | - |

### API Usage Examples

```bash
# Get pathway information
curl "https://rest.kegg.jp/get/hsa04010"

# List all human pathways
curl "https://rest.kegg.jp/list/pathway/hsa"

# Convert UniProt IDs to KEGG genes
curl "https://rest.kegg.jp/conv/hsa/uniprot:P00533"

# Get drug information
curl "https://rest.kegg.jp/get/D00001"
```

### Update Frequency

- **Weekly updates** (for subscription FTP)
- REST API reflects current database state

### Processing Recommendations

**For bulk access without subscription:**
1. Use REST API with appropriate delays between requests
2. Cache results locally to minimize repeated queries
3. Consider academic subscription for large-scale projects

```python
import time
import requests

def kegg_batch_get(ids, delay=0.5):
    """Batch retrieve KEGG entries with rate limiting"""
    results = []
    # KEGG allows max 10 IDs per request
    for i in range(0, len(ids), 10):
        batch = ids[i:i+10]
        url = f"https://rest.kegg.jp/get/{'+'.join(batch)}"
        response = requests.get(url)
        results.append(response.text)
        time.sleep(delay)  # Rate limiting
    return results
```

**Important Notes:**
- No official rate limits documented, but be respectful
- Image and KGML retrieval limited to single entries
- For comprehensive pathway analysis, academic subscription recommended

---

## 3. WikiPathways

**Website:** https://www.wikipathways.org/
**Data Archive:** https://data.wikipathways.org/
**License:** Creative Commons Attribution 4.0

### Monthly Release Schedule

- Releases on the **10th of each month**
- 12-month rolling archive maintained on data.wikipathways.org
- Zenodo archive for DOI-citable long-term storage

### Download URLs

| Format | URL | Description |
|--------|-----|-------------|
| GPML | https://data.wikipathways.org/current/gpml/ | Pathway markup language (XML-based) |
| GMT | https://data.wikipathways.org/current/gmt/ | Gene Matrix Transposed (gene sets) |
| SVG | https://data.wikipathways.org/current/svg/ | Scalable vector graphics |
| RDF | Via Zenodo | Semantic web format |

### Current Release Files (January 2026)

**GPML Files (Selected Organisms):**

| Organism | File | Size |
|----------|------|------|
| Homo sapiens | wikipathways-20260110-gpml-Homo_sapiens.zip | 11 MB |
| Mus musculus | wikipathways-20260110-gpml-Mus_musculus.zip | 1.9 MB |
| Bos taurus | wikipathways-20260110-gpml-Bos_taurus.zip | 2.5 MB |
| Rattus norvegicus | wikipathways-20260110-gpml-Rattus_norvegicus.zip | 928 KB |
| Saccharomyces cerevisiae | wikipathways-20260110-gpml-Saccharomyces_cerevisiae.zip | 484 KB |
| Danio rerio | wikipathways-20260110-gpml-Danio_rerio.zip | 412 KB |
| Arabidopsis thaliana | wikipathways-20260110-gpml-Arabidopsis_thaliana.zip | ~300 KB |

**GMT Files (Selected Organisms):**

| Organism | File | Size |
|----------|------|------|
| Homo sapiens | wikipathways-20260110-gmt-Homo_sapiens.gmt | 336 KB |
| Mus musculus | wikipathways-20260110-gmt-Mus_musculus.gmt | 92 KB |
| Bos taurus | wikipathways-20260110-gmt-Bos_taurus.gmt | 92 KB |
| Rattus norvegicus | wikipathways-20260110-gmt-Rattus_norvegicus.gmt | 48 KB |

### Programmatic Access

```bash
# Download current human GPML
wget https://data.wikipathways.org/current/gpml/wikipathways-*-gpml-Homo_sapiens.zip

# Download current human GMT
wget https://data.wikipathways.org/current/gmt/wikipathways-*-gmt-Homo_sapiens.gmt

# Download all GMT files
wget -r -np -nd -A "*.gmt" https://data.wikipathways.org/current/gmt/
```

**R Package (rWikiPathways):**

```R
library(rWikiPathways)

# Download latest human GMT
downloadPathwayArchive(organism="Homo sapiens", format="gmt")

# Download specific release
downloadPathwayArchive(organism="Homo sapiens", format="gmt", date="20260110")
```

### Processing Recommendations

- **GMT format** recommended for gene set enrichment analysis (GSEA, etc.)
- **GPML format** for detailed pathway visualization and editing
- Use Zenodo archives for reproducible research with DOIs
- Monthly updates mean relatively fresh data without excessive churn

---

## 4. UniProt

**Website:** https://www.uniprot.org/
**FTP:** https://ftp.uniprot.org/pub/databases/uniprot/
**License:** Creative Commons Attribution 4.0

### Knowledge Base Files

**Base URL:** https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/

#### Swiss-Prot (Reviewed, High-Quality)

| File | Size | Description |
|------|------|-------------|
| uniprot_sprot.xml.gz | 883 MB | Full XML format with all annotations |
| uniprot_sprot.fasta.gz | 89 MB | FASTA sequences only |
| uniprot_sprot.dat.gz | 655 MB | Flat file format |
| uniprot_sprot_varsplic.fasta.gz | 8.2 MB | Variant splice isoforms |

#### TrEMBL (Unreviewed, Comprehensive)

| File | Size | Description |
|------|------|-------------|
| uniprot_trembl.xml.gz | 171 GB | Full XML format |
| uniprot_trembl.fasta.gz | 48 GB | FASTA sequences only |
| uniprot_trembl.dat.gz | 139 GB | Flat file format |

### ID Mapping Files

**Base URL:** https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/

| File | Size | Description |
|------|------|-------------|
| idmapping.dat.gz | 18 GB | Complete ID mapping database |
| idmapping_selected.tab.gz | 9 GB | Curated subset of common mappings |
| by_organism/ | Varies | Species-specific mapping files |

### Update Frequency

- **Every 8 weeks** (approximately bi-monthly)
- Release notes at https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt

### Download Examples

```bash
# Download Swiss-Prot (recommended starting point)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz

# Download FASTA sequences
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# Download ID mappings for human only
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz

# Using rsync for efficient updates
rsync -avz rsync://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz .
```

### Processing Recommendations

**For most use cases:**
1. Start with **Swiss-Prot** (reviewed entries only) - more manageable size
2. Use **organism-specific ID mapping files** to reduce processing time
3. Consider **FASTA format** for sequence-only applications
4. Use **XML format** when full annotation details needed

```python
import gzip
from Bio import SeqIO

# Parse Swiss-Prot FASTA
with gzip.open('uniprot_sprot.fasta.gz', 'rt') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        # Process each protein
        print(record.id, len(record.seq))
```

**Storage Considerations:**
- Swiss-Prot: ~5 GB uncompressed (XML)
- TrEMBL: ~1 TB uncompressed (XML) - requires significant storage
- ID mappings: ~100 GB uncompressed

---

## 5. STRING

**Website:** https://string-db.org/
**Download Page:** https://string-db.org/cgi/download
**License:** Creative Commons Attribution 4.0

### Current Version: v12.0

### Protein Interaction Files

| File | Size | Description |
|------|------|-------------|
| protein.links.v12.0.txt.gz | 128.7 GB | Combined score protein links |
| protein.links.detailed.v12.0.txt.gz | 189.6 GB | Individual evidence channel scores |
| protein.links.full.v12.0.txt.gz | 199.6 GB | All scores including transferred |

### Physical Interaction Data

| File | Size | Description |
|------|------|-------------|
| protein.physical.links.v12.0.txt.gz | 11.1 GB | Physical binding interactions only |
| protein.physical.links.detailed.v12.0.txt.gz | 13.8 GB | Detailed physical scores |
| protein.physical.links.full.v12.0.txt.gz | 14.5 GB | Full physical interaction data |

### Protein Information Files

| File | Size | Description |
|------|------|-------------|
| protein.info.v12.0.txt.gz | 1.2 GB | Protein names and descriptions |
| protein.sequences.v12.0.fa.gz | 12.1 GB | Protein sequences (FASTA) |
| protein.aliases.v12.0.txt.gz | 3.2 GB | Alternative protein names/IDs |

### Species & Orthology Data

| File | Size | Description |
|------|------|-------------|
| species.v12.0.txt | 999 KB | List of organisms |
| species.tree.v12.0.txt | 92.8 MB | Taxonomic tree |
| COG.mappings.v12.0.txt.gz | 720 MB | Clusters of Orthologous Groups |

### Species-Specific Downloads

STRING provides per-species files for more manageable downloads:

```
https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
```

Where `9606` is the NCBI taxonomy ID for human.

### Download Examples

```bash
# Download human protein interactions only
wget https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz

# Download human protein info
wget https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz

# Download species list first to identify taxonomy IDs
wget https://stringdb-downloads.org/download/species.v12.0.txt
```

### Update Frequency

- **Major releases every 1-2 years**
- Current version: 12.0 (released 2023)

### Processing Recommendations

**For most applications:**
1. Download **species-specific files** rather than full database
2. Use **protein.links** for functional associations
3. Use **protein.physical.links** for direct binding evidence
4. Filter by **combined_score** threshold (typically >400 for medium confidence, >700 for high)

```python
import pandas as pd
import gzip

# Load human protein links
with gzip.open('9606.protein.links.v12.0.txt.gz', 'rt') as f:
    df = pd.read_csv(f, sep=' ')

# Filter for high confidence interactions
high_conf = df[df['combined_score'] >= 700]
print(f"High confidence interactions: {len(high_conf)}")
```

**Storage Note:** Full protein.links file (129 GB compressed) expands to ~1 TB uncompressed.

---

## 6. Gene Ontology

**Website:** http://geneontology.org/
**Download Page:** http://geneontology.org/docs/download-ontology/
**License:** Creative Commons Attribution 4.0

### Ontology Files

**Base URL:** https://purl.obolibrary.org/obo/go/

| File | URL | Description |
|------|-----|-------------|
| go-basic.obo | https://purl.obolibrary.org/obo/go/go-basic.obo | Filtered, acyclic version (recommended) |
| go.obo | https://purl.obolibrary.org/obo/go.obo | Full ontology in OBO format |
| go.owl | https://purl.obolibrary.org/obo/go.owl | Full ontology in OWL format |
| go.json | https://purl.obolibrary.org/obo/go.json | JSON format |

### Gene Annotation Files (GAF Format)

**Base URL:** https://current.geneontology.org/annotations/

| Organism | File | Annotations |
|----------|------|-------------|
| Human | goa_human.gaf.gz | 972,445 |
| Mouse | mgi.gaf.gz | 744,389 |
| Rat | rgd.gaf.gz | 590,055 |
| Arabidopsis | tair.gaf.gz | 227,211 |
| Drosophila | fb.gaf.gz | 186,975 |
| Yeast | sgd.gaf.gz | 169,668 |
| C. elegans | wb.gaf.gz | 125,157 |
| Zebrafish | zfin.gaf.gz | 250,491 |

### Additional Annotation Sources

| Source | URL | Description |
|--------|-----|-------------|
| UniProt GAFs | https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/ | Per-proteome annotations (~20,000 species) |
| NCBI RefSeq | Via NCBI Datasets | Reference genome annotations |
| Zenodo Archive | https://zenodo.org/record/1205166 | DOI-versioned monthly releases |

### Update Frequency

- **Monthly releases**
- Zenodo archive from August 2018 onward

### Download Examples

```bash
# Download ontology
wget https://purl.obolibrary.org/obo/go/go-basic.obo

# Download human annotations
wget https://current.geneontology.org/annotations/goa_human.gaf.gz

# Download all model organism annotations
for org in goa_human mgi rgd sgd fb wb zfin tair; do
    wget "https://current.geneontology.org/annotations/${org}.gaf.gz"
done
```

### Processing Recommendations

**Required files for most analyses:**
1. **go-basic.obo** - Ontology structure
2. **Species-specific GAF** - Gene annotations

```python
from goatools.obo_parser import GODag
from goatools.associations import read_gaf

# Load ontology
go_dag = GODag('go-basic.obo')

# Load human annotations
associations = read_gaf('goa_human.gaf.gz')
```

**Best Practices:**
- Use **go-basic.obo** for annotation tools (filtered, acyclic)
- Use **go.owl** for semantic web applications
- Always pair annotation files with matching ontology release

---

## 7. NCBI Gene

**Website:** https://www.ncbi.nlm.nih.gov/gene
**FTP:** https://ftp.ncbi.nlm.nih.gov/gene/DATA/
**License:** Public Domain

### Core Data Files

| File | Size | Description | Updated |
|------|------|-------------|---------|
| gene_info.gz | 1.3 GB | Primary gene annotation data | Daily |
| gene2refseq.gz | 2.0 GB | Gene-to-RefSeq sequence mappings | Daily |
| gene2go.gz | 1.2 GB | Gene Ontology associations | Daily |
| gene2accession.gz | 3.9 GB | Gene-to-GenBank accession mappings | Daily |
| gene2pubmed.gz | 237 MB | Gene-to-PubMed literature links | Daily |
| gene2ensembl.gz | 278 MB | Gene-to-Ensembl cross-references | Daily |
| gene_history.gz | 149 MB | Historical gene name/ID changes | Daily |
| gene_orthologs.gz | 107 MB | Orthologous gene relationships | Daily |
| gene_summary.gz | 21 MB | Concise gene descriptions | Daily |

### Species-Specific Files

Available in subdirectory structure:
```
https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/
https://ftp.ncbi.nlm.nih.gov/gene/DATA/ASN_BINARY/
```

### Update Frequency

- **Daily updates**
- Most comprehensive and current gene information source

### Download Examples

```bash
# Download core files
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz

# Filter for human only (taxonomy ID 9606)
zcat gene_info.gz | awk -F'\t' '$1==9606' > human_gene_info.txt
```

### Processing Recommendations

**Key fields in gene_info:**
- Column 1: tax_id
- Column 2: GeneID
- Column 3: Symbol
- Column 5: Synonyms
- Column 9: description
- Column 10: type_of_gene

```python
import pandas as pd
import gzip

# Load human genes
with gzip.open('gene_info.gz', 'rt') as f:
    # Skip comment lines
    df = pd.read_csv(f, sep='\t', comment='#',
                     names=['tax_id', 'GeneID', 'Symbol', 'LocusTag',
                            'Synonyms', 'dbXrefs', 'chromosome', 'map_location',
                            'description', 'type_of_gene', 'Symbol_from_nomenclature',
                            'Full_name_from_nomenclature', 'Nomenclature_status',
                            'Other_designations', 'Modification_date', 'Feature_type'])

    human_genes = df[df['tax_id'] == 9606]
    print(f"Human genes: {len(human_genes)}")
```

**Best Practices:**
- Filter by taxonomy ID early to reduce memory usage
- Use gene2ensembl for cross-referencing with Ensembl
- gene_history important for tracking deprecated/merged gene IDs

---

## 8. Ensembl

**Website:** https://www.ensembl.org/
**FTP:** https://ftp.ensembl.org/pub/
**License:** Open access (varies by data type)

### Current Release: 115 (September 2025)

### FTP Directory Structure

```
https://ftp.ensembl.org/pub/release-115/
├── fasta/           # Sequence files
│   └── homo_sapiens/
│       ├── cdna/    # cDNA sequences
│       ├── cds/     # Coding sequences
│       ├── dna/     # Genomic DNA
│       ├── ncrna/   # Non-coding RNA
│       └── pep/     # Protein sequences
├── gtf/             # Gene annotations (GTF format)
├── gff3/            # Gene annotations (GFF3 format)
├── mysql/           # Database dumps
├── variation/       # Variant data (VCF, GVF)
└── regulation/      # Regulatory features
```

### Human GTF/GFF3 Files

**Base URL:** https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/

| File | Size | Description |
|------|------|-------------|
| Homo_sapiens.GRCh38.115.gtf.gz | 100 MB | Complete GTF annotation |
| Homo_sapiens.GRCh38.115.chr.gtf.gz | 99 MB | Chromosomes only |
| Homo_sapiens.GRCh38.115.chr_patch_hapl_scaff.gtf.gz | 104 MB | Including patches/haplotypes |
| Homo_sapiens.GRCh38.115.abinitio.gtf.gz | 3.2 MB | Ab initio predictions |

### Human FASTA Files

| Directory | Contents |
|-----------|----------|
| dna/ | Genomic DNA (primary assembly ~900 MB compressed) |
| cdna/ | All transcript sequences |
| cds/ | Coding sequences only |
| pep/ | Protein sequences |
| ncrna/ | Non-coding RNA sequences |

### BioMart Bulk Exports

**URL:** https://www.ensembl.org/biomart/martview

BioMart provides customized bulk exports with:
- User-defined attribute selection
- Multiple filter options
- Cross-database queries (Ensembl, UniProt, etc.)
- XML query format for programmatic access

### Update Frequency

- **3-4 releases per year** (approximately quarterly)
- Each release includes updated genome assemblies and annotations

### Download Methods

```bash
# Direct wget
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

# rsync for efficient updates
rsync -avz rsync://ftp.ensembl.org/ensembl/pub/release-115/gtf/homo_sapiens/ ./ensembl_gtf/

# Globus for large transfers
# Use endpoint: Ensembl Public Data
```

### BioMart Programmatic Access

```python
from biomart import BiomartServer

server = BiomartServer("http://www.ensembl.org/biomart")
ensembl = server.datasets['hsapiens_gene_ensembl']

# Query gene information
response = ensembl.search({
    'filters': {'chromosome_name': ['1', '2', '3']},
    'attributes': ['ensembl_gene_id', 'external_gene_name', 'description']
})

for line in response.iter_lines():
    print(line.decode('utf-8'))
```

### Processing Recommendations

**For gene annotation:**
- Use **GTF format** for most downstream tools (RNA-seq, etc.)
- **GFF3** for more detailed feature hierarchies

**For sequences:**
- **Primary assembly** DNA recommended (excludes patches/haplotypes)
- Download only needed sequence types (cdna/cds/pep/ncrna)

```bash
# Download human reference for RNA-seq
wget https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

---

## 9. Processing Recommendations Summary

### Storage Requirements

| Database | Minimal Install | Full Install |
|----------|-----------------|--------------|
| Reactome | ~100 MB (text files) | ~3 GB (with Neo4j) |
| KEGG | Via API only | Subscription required |
| WikiPathways | ~15 MB (human GMT) | ~50 MB (all species) |
| UniProt (Swiss-Prot) | ~100 MB (FASTA) | ~5 GB (XML) |
| UniProt (TrEMBL) | ~50 GB (FASTA) | ~1 TB (XML) |
| STRING | ~15 GB (human only) | ~200 GB (all species) |
| Gene Ontology | ~50 MB | ~500 MB (with all annotations) |
| NCBI Gene | ~500 MB (human) | ~10 GB (all species) |
| Ensembl | ~200 MB (human GTF+FASTA) | ~50 GB (full release) |

### Recommended Download Priority

**For pathway analysis:**
1. WikiPathways GMT files (smallest, monthly updates)
2. Reactome mapping files
3. KEGG via API (as needed)

**For protein/gene mapping:**
1. UniProt Swiss-Prot (quality over quantity)
2. NCBI gene_info + gene2ensembl
3. UniProt ID mapping files

**For protein interactions:**
1. STRING species-specific files
2. Filter by confidence threshold

**For functional annotation:**
1. GO ontology (go-basic.obo)
2. Species-specific GAF files
3. Ensembl GTF for genomic context

### Update Schedule Recommendation

| Database | Check Frequency | Rationale |
|----------|-----------------|-----------|
| NCBI Gene | Weekly | Daily updates |
| Gene Ontology | Monthly | Monthly releases |
| WikiPathways | Monthly | Monthly releases |
| UniProt | Bi-monthly | 8-week release cycle |
| Reactome | Quarterly | Quarterly releases |
| Ensembl | Quarterly | 3-4 releases/year |
| STRING | Yearly | Major releases every 1-2 years |

### Automation Script Template

```bash
#!/bin/bash
# bulk_download.sh - Download core biological databases

DATA_DIR="/data/databases"
mkdir -p $DATA_DIR/{reactome,wikipathways,uniprot,ncbi,go,ensembl,string}

# Reactome
echo "Downloading Reactome..."
cd $DATA_DIR/reactome
wget -N https://reactome.org/download/current/UniProt2Reactome.txt
wget -N https://reactome.org/download/current/NCBI2Reactome.txt
wget -N https://reactome.org/download/current/ReactomePathways.txt

# WikiPathways
echo "Downloading WikiPathways..."
cd $DATA_DIR/wikipathways
wget -N https://data.wikipathways.org/current/gmt/wikipathways-*-gmt-Homo_sapiens.gmt

# UniProt (Swiss-Prot only)
echo "Downloading UniProt Swiss-Prot..."
cd $DATA_DIR/uniprot
wget -N https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# NCBI Gene
echo "Downloading NCBI Gene..."
cd $DATA_DIR/ncbi
wget -N https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
wget -N https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
wget -N https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz

# Gene Ontology
echo "Downloading Gene Ontology..."
cd $DATA_DIR/go
wget -N https://purl.obolibrary.org/obo/go/go-basic.obo
wget -N https://current.geneontology.org/annotations/goa_human.gaf.gz

# Ensembl (human)
echo "Downloading Ensembl human..."
cd $DATA_DIR/ensembl
wget -N https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.*.gtf.gz

echo "Download complete!"
```

---

## References

- Reactome: https://reactome.org/download-data
- KEGG: https://www.kegg.jp/kegg/rest/keggapi.html
- WikiPathways: https://data.wikipathways.org/
- UniProt: https://ftp.uniprot.org/pub/databases/uniprot/
- STRING: https://string-db.org/cgi/download
- Gene Ontology: http://geneontology.org/docs/download-ontology/
- NCBI Gene: https://ftp.ncbi.nlm.nih.gov/gene/DATA/
- Ensembl: https://ftp.ensembl.org/pub/
