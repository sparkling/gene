---
id: schema-biogrid
title: "BioGRID Schema Documentation"
type: schema
parent: _index.md
last_updated: 2026-01-24
status: final
tags: [schema, database, protein-interactions, genetic-interactions]
---

# BioGRID - Schema Documentation

## TL;DR

BioGRID provides protein-protein and genetic interaction data in multiple tab-delimited formats. The primary format is TAB 2.0/3.0, with PSI-MITAB available for standard interchange.

## TAB 2.0 Format (26 Columns)

### Column Definitions

| Column | Name | Description |
|--------|------|-------------|
| 1 | BioGRID Interaction ID | Unique interaction identifier |
| 2 | Entrez Gene ID Interactor A | NCBI Gene ID for interactor A |
| 3 | Entrez Gene ID Interactor B | NCBI Gene ID for interactor B |
| 4 | BioGRID ID Interactor A | BioGRID internal ID for A |
| 5 | BioGRID ID Interactor B | BioGRID internal ID for B |
| 6 | Systematic Name Interactor A | Systematic gene name for A |
| 7 | Systematic Name Interactor B | Systematic gene name for B |
| 8 | Official Symbol Interactor A | HGNC/official symbol for A |
| 9 | Official Symbol Interactor B | HGNC/official symbol for B |
| 10 | Synonyms Interactor A | Alias names (pipe-separated) |
| 11 | Synonyms Interactor B | Alias names (pipe-separated) |
| 12 | Experimental System | Detection method |
| 13 | Experimental System Type | physical or genetic |
| 14 | Author | First author surname |
| 15 | Publication Source | PubMed ID |
| 16 | Organism ID Interactor A | NCBI Taxonomy ID for A |
| 17 | Organism ID Interactor B | NCBI Taxonomy ID for B |
| 18 | Throughput | high or low |
| 19 | Score | Confidence score (if available) |
| 20 | Modification | PTM information |
| 21 | Qualifications | Additional details |
| 22 | Tags | Custom tags |
| 23 | Source Database | Source of curation |
| 24 | SWISS-PROT Accession Interactor A | UniProt AC for A |
| 25 | TREMBL Accession Interactor A | TrEMBL AC for A |
| 26 | SWISS-PROT Accession Interactor B | UniProt AC for B |

## TAB 3.0 Format (35 Columns)

Additional columns in TAB 3.0:

| Column | Name | Description |
|--------|------|-------------|
| 27 | TREMBL Accession Interactor B | TrEMBL AC for B |
| 28 | REFSEQ Accession Interactor A | RefSeq AC for A |
| 29 | REFSEQ Accession Interactor B | RefSeq AC for B |
| 30 | Ontology Term IDs | GO terms (pipe-separated) |
| 31 | Ontology Term Names | GO term names |
| 32 | Ontology Term Categories | BP, MF, CC |
| 33 | Ontology Term Qualifier IDs | Qualifier IDs |
| 34 | Ontology Term Qualifier Names | Qualifier names |
| 35 | Ontology Term Types | Term types |

## Example TAB 2.0 Record

```
103	7157	3845	112315	107485	-	-	TP53	KRAS	p53|LFS1	K-RAS|RASK2	Two-hybrid	physical	Vetter	1535557	9606	9606	Low Throughput	-	-	-	-	BioGRID	P04637	-	P01116
```

## PSI-MITAB 2.5 Format

Standard 15-column format:

| Column | Field | Example |
|--------|-------|---------|
| 1 | ID Interactor A | entrez gene:7157 |
| 2 | ID Interactor B | entrez gene:3845 |
| 3 | Alt ID A | biogrid:112315 |
| 4 | Alt ID B | biogrid:107485 |
| 5 | Alias A | entrez gene:TP53(gene name) |
| 6 | Alias B | entrez gene:KRAS(gene name) |
| 7 | Detection Method | psi-mi:"MI:0018"(two hybrid) |
| 8 | First Author | Vetter IR (1997) |
| 9 | Publication | pubmed:1535557 |
| 10 | Taxon A | taxid:9606(human) |
| 11 | Taxon B | taxid:9606(human) |
| 12 | Interaction Type | psi-mi:"MI:0407"(direct interaction) |
| 13 | Source DB | psi-mi:"MI:0463"(biogrid) |
| 14 | Interaction ID | biogrid:103 |
| 15 | Confidence | - |

## Experimental Systems

### Physical Interactions

| System | MI ID | Description |
|--------|-------|-------------|
| Two-hybrid | MI:0018 | Yeast two-hybrid |
| Affinity Capture-MS | MI:0004 | Mass spec pull-down |
| Affinity Capture-Western | MI:0004 | Western blot validation |
| Co-fractionation | MI:0027 | Co-purification |
| Co-crystal Structure | MI:0114 | X-ray crystallography |
| Reconstituted Complex | MI:0045 | In vitro reconstitution |
| PCA | MI:0090 | Protein complementation |
| FRET | MI:0055 | Fluorescence resonance |
| Proximity Label-MS | MI:1313 | BioID, APEX |
| Co-localization | MI:0403 | Microscopy co-localization |

### Genetic Interactions

| System | Description |
|--------|-------------|
| Synthetic Lethality | Lethal when both genes disrupted |
| Synthetic Growth Defect | Growth impairment |
| Synthetic Rescue | Suppressor interaction |
| Dosage Lethality | Gene dosage effect |
| Dosage Rescue | Dosage suppression |
| Phenotypic Enhancement | Phenotype worsening |
| Phenotypic Suppression | Phenotype improvement |
| Negative Genetic | Aggravating |
| Positive Genetic | Alleviating |

## Modification Types

| Modification | Description |
|--------------|-------------|
| Phosphorylation | Phosphorylated residue |
| Ubiquitination | Ubiquitin conjugation |
| Acetylation | Acetyl group addition |
| Methylation | Methyl group addition |
| Sumoylation | SUMO conjugation |

## JSON API Response

```json
{
  "interactions": [
    {
      "INTERACTION_ID": 103,
      "ENTREZ_GENE_A": "7157",
      "ENTREZ_GENE_B": "3845",
      "BIOGRID_A": 112315,
      "BIOGRID_B": 107485,
      "OFFICIAL_SYMBOL_A": "TP53",
      "OFFICIAL_SYMBOL_B": "KRAS",
      "EXPERIMENTAL_SYSTEM": "Two-hybrid",
      "EXPERIMENTAL_SYSTEM_TYPE": "physical",
      "PUBMED_ID": 1535557,
      "ORGANISM_A": 9606,
      "ORGANISM_B": 9606,
      "THROUGHPUT": "Low Throughput",
      "QUANTITATION": null,
      "MODIFICATION": null,
      "SOURCE_DATABASE": "BioGRID"
    }
  ]
}
```

## Organism Taxonomy IDs

| Organism | Tax ID |
|----------|--------|
| Homo sapiens | 9606 |
| Mus musculus | 10090 |
| Saccharomyces cerevisiae | 559292 |
| Drosophila melanogaster | 7227 |
| Caenorhabditis elegans | 6239 |
| Arabidopsis thaliana | 3702 |
| Schizosaccharomyces pombe | 284812 |
| Escherichia coli | 83333 |

## See Also

- [Download Documentation](./download.md)
- [BioGRID Downloads](https://downloads.thebiogrid.org/)
