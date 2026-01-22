# UniProt ID Mapping Schema

**Document ID:** UNIPROT-IDMAPPING-SCHEMA
**Status:** Final
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [xrefs.md](../integration/xrefs.md)

---

## Overview

UniProt ID Mapping is the master cross-reference service linking 286 databases to UniProt protein entries. It provides the most comprehensive protein-centric ID mapping available, with pre-computed files updated every 8 weeks synchronized with UniProt releases.

### Service Endpoints

| Resource | URL |
|----------|-----|
| **Web Interface** | https://www.uniprot.org/id-mapping |
| **REST API** | https://rest.uniprot.org/idmapping/ |
| **FTP Downloads** | ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ |
| **Field Configuration** | https://rest.uniprot.org/configure/idmapping/fields |

---

## File Formats

### idmapping.dat (~15 GB compressed)

Three-column tab-delimited file containing all mappings.

| Column | Name | Description | Example |
|--------|------|-------------|---------|
| 1 | UniProtKB-AC | UniProt accession | P04637 |
| 2 | ID_type | External database name | GeneID |
| 3 | ID | External identifier | 7157 |

**Sample Records:**
```
P04637	UniProtKB-ID	P53_HUMAN
P04637	GeneID	7157
P04637	RefSeq	NP_000537.3
P04637	RefSeq_NT	NM_000546.6
P04637	PDB	1AIE
P04637	PDB	1C26
P04637	GO	GO:0000785
P04637	GO	GO:0001046
P04637	Ensembl	ENST00000269305
P04637	HGNC	HGNC:11998
```

### idmapping_selected.tab (~3 GB compressed)

Pre-computed 22-column table with most frequently requested mappings.

| Col | Field Name | Description | Format | Example |
|-----|------------|-------------|--------|---------|
| 1 | UniProtKB-AC | UniProt accession | [A-Z][0-9][A-Z0-9]{3}[0-9] | P04637 |
| 2 | UniProtKB-ID | UniProt entry name | *_SPECIES | P53_HUMAN |
| 3 | GeneID | NCBI Entrez Gene | Numeric | 7157 |
| 4 | RefSeq | RefSeq protein IDs | NP_*; separated | NP_000537.3 |
| 5 | GI | GenInfo Identifier | Numeric; deprecated | 120407068 |
| 6 | PDB | Protein Data Bank | 4-char; separated | 1AIE;1C26;1DT7 |
| 7 | GO | Gene Ontology terms | GO:[0-9]{7}; separated | GO:0000785;GO:0001046 |
| 8 | UniRef100 | UniRef100 cluster | UniRef100_* | UniRef100_P04637 |
| 9 | UniRef90 | UniRef90 cluster | UniRef90_* | UniRef90_P04637 |
| 10 | UniRef50 | UniRef50 cluster | UniRef50_* | UniRef50_P04637 |
| 11 | UniParc | UniParc ID | UPI[0-9A-F]{10} | UPI000002ED67 |
| 12 | PIR | PIR accession | Letter+digits | A25224 |
| 13 | NCBI-taxon | Taxonomy ID | Numeric | 9606 |
| 14 | MIM | OMIM ID | 6 digits | 191170 |
| 15 | UniGene | UniGene cluster | Hs.[0-9]+ | Hs.408312 |
| 16 | PubMed | PubMed IDs | Numeric; separated | 11433014;11488916 |
| 17 | EMBL | EMBL accession | Alphanumeric; separated | AB082923;AF052180 |
| 18 | EMBL-CDS | EMBL CDS IDs | *.[0-9]; separated | BAC16799.1;AAC12971.1 |
| 19 | Ensembl | Ensembl transcript | ENST[0-9]{11} | ENST00000269305 |
| 20 | Ensembl_TRS | Ensembl transcript | ENST[0-9]{11} | ENST00000269305 |
| 21 | Ensembl_PRO | Ensembl protein | ENSP[0-9]{11} | ENSP00000269305 |
| 22 | Additional_PubMed | Extra PubMed refs | Numeric; separated | Additional references |

**Notes:**
- Multiple values within a field are semicolon-separated
- Empty fields are represented by empty strings
- File is sorted by UniProtKB-AC

**Sample Row (TP53):**
```
P04637	P53_HUMAN	7157	NP_000537.3	120407068	1AIE;1C26;1DT7;1GZH	GO:0000785;GO:0001046;GO:0001085	UniRef100_P04637	UniRef90_P04637	UniRef50_P04637	UPI000002ED67	A25224	9606	191170	Hs.408312	11433014;11488916;11522829	AB082923;AF052180;AF307851	BAC16799.1;AAC12971.1;AAG34762.1	ENST00000269305	ENST00000269305	ENSP00000269305
```

---

## Complete Database List (286 Databases)

### Sequence Databases (8)

| Database | ID Type | Description |
|----------|---------|-------------|
| EMBL | EMBL/GenBank/DDBJ | Nucleotide sequence |
| GenBank | GenBank | NCBI nucleotide |
| DDBJ | DDBJ | DNA Data Bank Japan |
| RefSeq | RefSeq Protein | NCBI RefSeq protein |
| RefSeq_NT | RefSeq Nucleotide | NCBI RefSeq nucleotide |
| CCDS | CCDS | Consensus CDS project |
| PIR | PIR | Protein Information Resource |
| GI | GI number | GenInfo Identifier (deprecated) |

### Structure Databases (5)

| Database | ID Type | Description |
|----------|---------|-------------|
| PDB | PDB | Protein Data Bank |
| AlphaFoldDB | AF | AlphaFold predictions |
| SMR | SMR | Swiss-Model Repository |
| ModBase | ModBase | Comparative models |
| BMRB | BMRB | NMR data |

### Genome Annotation (12)

| Database | ID Type | Description |
|----------|---------|-------------|
| Ensembl | ENSG/ENST/ENSP | Ensembl gene/transcript/protein |
| Ensembl_TRS | ENST | Ensembl transcript |
| Ensembl_PRO | ENSP | Ensembl protein |
| GeneID | Numeric | NCBI Entrez Gene |
| KEGG | K/ko numbers | KEGG gene/ortholog |
| UCSC | Various | UCSC Genome Browser |
| PATRIC | PATRIC ID | Bacterial genome |
| WBParaSite | WBGene | Parasite genomes |
| MANE-Select | Various | Matched annotation |
| Ensembl_Genomes | Various | Non-vertebrate Ensembl |
| GeneDB | Various | Pathogen genomes |
| VEuPathDB | Various | Eukaryotic pathogens |

### Protein Interactions (10)

| Database | ID Type | Description |
|----------|---------|-------------|
| BioGRID | Numeric | Interaction database |
| ComplexPortal | CPX | Protein complexes |
| DIP | DIP | Database of Interacting Proteins |
| STRING | STRING | Protein networks |
| IntAct | EBI | Molecular interactions |
| MINT | MINT | Molecular Interaction Database |
| CORUM | CORUM | Mammalian complexes |
| ELM | ELM | Linear motifs |
| FunCoup | FunCoup | Functional coupling |
| mentha | mentha | Protein interactions |

### Chemistry Databases (8)

| Database | ID Type | Description |
|----------|---------|-------------|
| ChEMBL | CHEMBL | Bioactivity data |
| DrugBank | DB | Drug information |
| BindingDB | Numeric | Binding affinity |
| SwissLipids | SLM | Lipid database |
| GuidetoPHARMACOLOGY | Numeric | Pharmacology |
| DrugCentral | DC | Drug information |
| ZINC | ZINC | Commercial compounds |
| PubChem | CID | Chemical structures |

### Family & Domain Databases (15)

| Database | ID Type | Description |
|----------|---------|-------------|
| InterPro | IPR | Protein signatures |
| Pfam | PF | Protein families |
| SMART | SM | Domain architecture |
| PROSITE | PS | Patterns/profiles |
| Gene3D | G3D | Structural domains |
| HAMAP | MF | Microbial annotation |
| CDD | cd | Conserved domains |
| Superfamily | SSF | Structural classification |
| DisProt | DP | Disordered proteins |
| MobiDB | Various | Mobility annotation |
| PRINTS | PR | Fingerprints |
| PIRSF | PIRSF | Protein families |
| PANTHER | PTHR | Classification |
| TIGRFAMs | TIGR | Microbial families |
| NCBIfam | NF | NCBI families |

### Organism-Specific Databases (38)

| Database | Organism | Description |
|----------|----------|-------------|
| HGNC | Human | Gene nomenclature |
| MGI | Mouse | Mouse Genome Informatics |
| RGD | Rat | Rat Genome Database |
| FlyBase | Drosophila | Fly genome |
| WormBase | C. elegans | Worm genome |
| SGD | S. cerevisiae | Yeast genome |
| PomBase | S. pombe | Fission yeast |
| TAIR | Arabidopsis | Plant genome |
| ZFIN | Zebrafish | Fish genome |
| Xenbase | Xenopus | Frog genome |
| Araport | Arabidopsis | Plant portal |
| CTD | Various | Comparative toxicogenomics |
| neXtProt | Human | Human proteome |
| HPA | Human | Protein Atlas |
| dictyBase | Dictyostelium | Slime mold |
| EcoGene | E. coli | Bacterial genome |
| CGD | Candida | Fungal genome |
| VectorBase | Vectors | Disease vectors |
| EuPathDB | Pathogens | Eukaryotic pathogens |
| GeneCards | Human | Gene database |
| MalaCards | Human | Disease database |
| (+ 17 more organism databases) | | |

### PTM Databases (10)

| Database | ID Type | Description |
|----------|---------|-------------|
| PhosphoSitePlus | PSP | Phosphorylation |
| CarbonylDB | CDB | Carbonylation |
| DEPOD | DEPOD | Phosphatases |
| iPTMnet | iPTMnet | PTM network |
| SwissPalm | SwissPalm | Palmitoylation |
| GlyConnect | GlyConnect | Glycosylation |
| GlyCosmos | GlyCosmos | Glycoscience |
| GlyGen | GlyGen | Glycoinformatics |
| MetOSite | MetOSite | Oxidation |
| PhosphoNET | Various | Phosphorylation |

### Enzyme & Pathway (15)

| Database | ID Type | Description |
|----------|---------|-------------|
| BioCyc | Various | Pathway/genome |
| BRENDA | EC | Enzyme database |
| Reactome | R-HSA | Pathways |
| SABIO-RK | SABIO | Reaction kinetics |
| UniPathway | UPA | Metabolic pathways |
| SignaLink | SL | Signaling |
| SIGNOR | SIGNOR | Signaling |
| PlantReactome | Various | Plant pathways |
| PathBank | PW | Pathways |
| CAZy | Various | Carbohydrate enzymes |
| ESTHER | ESTHER | Esterases |
| MoonDB | MoonDB | Moonlighting |
| PeroxiBase | PB | Peroxidases |
| REBASE | REBASE | Restriction enzymes |
| KEGG | K/map | Pathways |

### Proteomic Databases (10)

| Database | ID Type | Description |
|----------|---------|-------------|
| PRIDE | PXD | Proteomics data |
| PeptideAtlas | Various | Peptide atlas |
| PaxDb | Various | Protein abundance |
| ProteomicsDB | Various | MS data |
| MassIVE | MSV | Mass spectrometry |
| jPOST | JPST | Japanese proteome |
| CPTAC | Various | Cancer proteome |
| TopDownProteomics | Various | Top-down MS |
| Pumba | Various | MS data |
| ProMEX | Various | Plant MS |

### Phylogenomic (8)

| Database | ID Type | Description |
|----------|---------|-------------|
| eggNOG | ENOG | Ortholog groups |
| OMA | OMA | Ortholog matrix |
| OrthoDB | OrthoDB | Ortholog database |
| InParanoid | Various | Orthologs |
| GeneTree | ENSGT | Gene trees |
| TreeFam | TF | Gene families |
| PhylomeDB | Phy | Phylomes |
| PAN-GO | Various | Pan genome GO |

### Genetic Variation (6)

| Database | ID Type | Description |
|----------|---------|-------------|
| dbSNP | rs | SNP database |
| ClinGen | Various | Clinical genome |
| GenCC | GenCC | Gene curation |
| DisGeNET | Various | Gene-disease |
| DMDM | Various | Domain mutations |
| BioMuta | Various | Cancer mutations |

### Gene Ontology & Function (5)

| Database | ID Type | Description |
|----------|---------|-------------|
| GO | GO:[0-9]{7} | Gene Ontology |
| QuickGO | GO | GO annotations |
| AmiGO | GO | GO browser |
| Reactome | R-HSA | Pathway GO |
| UniPathway | UPA | Metabolic GO |

---

## REST API Specification

### Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/run` | POST | Submit mapping job |
| `/status/{jobId}` | GET | Check job status |
| `/results/{jobId}` | GET | Get results |
| `/stream/{jobId}` | GET | Stream large results |

### Submit Job Request

```bash
curl -X POST "https://rest.uniprot.org/idmapping/run" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "from=UniProtKB_AC-ID&to=GeneID&ids=P04637,P01308,P00533"
```

**Response:**
```json
{
  "jobId": "a1b2c3d4e5f6g7h8i9j0"
}
```

### Check Status

```bash
curl "https://rest.uniprot.org/idmapping/status/a1b2c3d4e5f6g7h8i9j0"
```

**Response:**
```json
{
  "jobStatus": "FINISHED"
}
```

### Get Results

```bash
curl "https://rest.uniprot.org/idmapping/results/a1b2c3d4e5f6g7h8i9j0"
```

**Response:**
```json
{
  "results": [
    {
      "from": "P04637",
      "to": {
        "primaryAccession": "7157",
        "uniProtkbId": "P53_HUMAN"
      }
    }
  ]
}
```

### Rate Limits

| Access Type | Rate |
|-------------|------|
| Anonymous | 1 request/second |
| With email header | 3 requests/second |
| Bulk download | Preferred method |

---

## Supported ID Types for Mapping

### From UniProt

```
UniProtKB
UniProtKB_AC-ID
UniProtKB/Swiss-Prot
UniParc
UniRef50
UniRef90
UniRef100
```

### To/From External

```
CCDS
EMBL/GenBank/DDBJ
EMBL/GenBank/DDBJ CDS
Ensembl
Ensembl Genomes
Ensembl Genomes Protein
Ensembl Genomes Transcript
Ensembl Protein
Ensembl Transcript
Gene Name
GeneID
GI number
KEGG
PATRIC
RefSeq Nucleotide
RefSeq Protein
UCSC
WBParaSite
WBParaSite Transcript/Protein
```

### Organism-Specific

```
ArachnoServer
Araport
CGD
ConoServer
dictyBase
EchoBASE
euHCVdb
FlyBase
GeneCards
GeneReviews
HGNC
LegioList
Leproma
MaizeGDB
MGI
MIM
neXtProt
OpenTargets
Orphanet
PharmGKB
PomBase
PseudoCAP
RGD
SGD
TubercuList
VEuPathDB
VGNC
WormBase
WormBase Protein
WormBase Transcript
Xenbase
ZFIN
```

---

## Download Files

### FTP Structure

```
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
├── idmapping.dat.gz                    # ~15 GB compressed
├── idmapping_selected.tab.gz           # ~3 GB compressed
├── README                              # Documentation
└── by_organism/
    ├── HUMAN_9606_idmapping.dat.gz     # Human only
    ├── MOUSE_10090_idmapping.dat.gz    # Mouse only
    ├── RAT_10116_idmapping.dat.gz      # Rat only
    └── ...                             # Other organisms
```

### File Sizes (Approximate)

| File | Compressed | Uncompressed |
|------|------------|--------------|
| idmapping.dat.gz | ~15 GB | ~150 GB |
| idmapping_selected.tab.gz | ~3 GB | ~30 GB |
| HUMAN_9606_idmapping.dat.gz | ~500 MB | ~5 GB |

---

## Integration Examples

### Python - Parse Selected Tab File

```python
import gzip
import csv

COLUMNS = [
    'UniProtKB_AC', 'UniProtKB_ID', 'GeneID', 'RefSeq', 'GI',
    'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc',
    'PIR', 'NCBI_taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL',
    'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional_PubMed'
]

def parse_idmapping_selected(filepath):
    with gzip.open(filepath, 'rt') as f:
        reader = csv.DictReader(f, fieldnames=COLUMNS, delimiter='\t')
        for row in reader:
            # Parse multi-value fields
            if row['GO']:
                row['GO'] = row['GO'].split(';')
            if row['PDB']:
                row['PDB'] = row['PDB'].split(';')
            if row['PubMed']:
                row['PubMed'] = row['PubMed'].split(';')
            yield row

# Usage
for entry in parse_idmapping_selected('idmapping_selected.tab.gz'):
    if entry['UniProtKB_AC'] == 'P04637':
        print(f"TP53: GeneID={entry['GeneID']}, GO={entry['GO'][:3]}")
```

### Build Cross-Reference Index

```python
import sqlite3
import gzip

def build_index(filepath, db_path):
    conn = sqlite3.connect(db_path)
    conn.execute('''
        CREATE TABLE IF NOT EXISTS id_mapping (
            uniprot_ac TEXT PRIMARY KEY,
            gene_id TEXT,
            ensembl TEXT,
            refseq TEXT,
            pdb TEXT,
            go_terms TEXT,
            hgnc TEXT
        )
    ''')

    conn.execute('CREATE INDEX IF NOT EXISTS idx_gene_id ON id_mapping(gene_id)')
    conn.execute('CREATE INDEX IF NOT EXISTS idx_ensembl ON id_mapping(ensembl)')

    with gzip.open(filepath, 'rt') as f:
        reader = csv.DictReader(f, fieldnames=COLUMNS, delimiter='\t')
        batch = []
        for row in reader:
            batch.append((
                row['UniProtKB_AC'],
                row['GeneID'],
                row['Ensembl'],
                row['RefSeq'],
                row['PDB'],
                row['GO'],
                None  # HGNC needs separate lookup
            ))
            if len(batch) >= 10000:
                conn.executemany(
                    'INSERT OR REPLACE INTO id_mapping VALUES (?,?,?,?,?,?,?)',
                    batch
                )
                batch = []
        if batch:
            conn.executemany(
                'INSERT OR REPLACE INTO id_mapping VALUES (?,?,?,?,?,?,?)',
                batch
            )
    conn.commit()
    conn.close()
```

---

## Statistics

| Metric | Value |
|--------|-------|
| Total databases | 286 |
| UniProtKB entries | ~250 million |
| Swiss-Prot (reviewed) | ~570,000 |
| TrEMBL (unreviewed) | ~250 million |
| Human proteins | ~20,000 (Swiss-Prot) |
| Update frequency | Every 8 weeks |
| License | CC BY 4.0 |

---

## References

- UniProt Consortium (2025) "UniProt: the Universal Protein Knowledgebase in 2025" Nucleic Acids Research
- FTP README: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README
- API Documentation: https://www.uniprot.org/help/api_idmapping
