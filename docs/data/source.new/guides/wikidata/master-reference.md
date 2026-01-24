---
id: guides-wikidata-master-reference
title: "Complete Wikidata Extraction Guide for Genetics/Health Knowledge Base"
type: guide
parent: _index.md
last_updated: 2026-01-23
status: active
tags: [wikidata, sparql, extraction, knowledge-graph, q-ids, properties, guide]
---

**Parent:** [Wikidata Guides](./_index.md)

# Complete Wikidata Extraction Guide for Genetics/Health Knowledge Base

**Last Updated:** January 2026
**Purpose:** Comprehensive reference for extracting ALL relevant biomedical data from Wikidata

---

## Table of Contents

1. [Complete Q-ID Reference](#1-complete-q-id-reference)
2. [Complete Property (P-code) Reference](#2-complete-property-p-code-reference)
3. [Dump Filtering Strategy](#3-dump-filtering-strategy)
4. [Master SPARQL Queries](#4-master-sparql-queries)
5. [Data Completeness Matrix](#5-data-completeness-matrix)
6. [Recommended Extraction Pipeline](#6-recommended-extraction-pipeline)
7. [Complete Wikidata Extractor Code](#7-complete-wikidata-extractor-code)

---

## 1. Complete Q-ID Reference

### 1.1 Drugs and Medications

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q12140** | medication | Pharmaceutical drug (primary class) | ~45,000 |
| **Q35456** | pharmaceutical drug | Alternative drug class | ~2,500 |
| **Q8386** | drug | General drug concept | ~1,000 |
| **Q422248** | monoclonal antibody | Biologic drugs | ~3,500 |
| **Q17143904** | small molecule drug | Small molecule therapeutics | ~8,000 |
| **Q27136782** | biological medicine | Biologics | ~2,000 |
| **Q18214535** | prodrug | Inactive precursors | ~800 |
| **Q79460** | vaccine | Immunizations | ~1,200 |
| **Q58624061** | generic drug | Generic medications | ~5,000 |
| **Q7866** | radioactive tracer | Diagnostic/therapeutic | ~300 |

### 1.2 Chemical Compounds

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q11173** | chemical compound | General chemical (very broad) | ~2,500,000 |
| **Q79529** | chemical substance | Substances | ~500,000 |
| **Q113145171** | type of chemical entity | Chemical entity types | ~10,000 |
| **Q11344** | alkaloid | Plant alkaloids | ~12,000 |
| **Q134219** | flavonoid | Plant flavonoids | ~8,000 |
| **Q190875** | terpenoid | Terpene derivatives | ~15,000 |
| **Q407595** | glycoside | Sugar derivatives | ~6,000 |
| **Q2725376** | polyphenol | Phenolic compounds | ~4,000 |
| **Q193893** | steroid | Steroid compounds | ~8,000 |
| **Q3314483** | coumarins | Coumarin compounds | ~2,000 |

### 1.3 Traditional Chinese Medicine (TCM)

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q891104** | traditional Chinese medicine | TCM as practice | ~100 |
| **Q51128287** | traditional Chinese medicine preparation | TCM formulas | ~5,000 |
| **Q18025426** | Chinese herbal medicine | Chinese herbs | ~3,000 |
| **Q54951330** | component of traditional Chinese medicine | TCM ingredients | ~2,000 |
| **Q899494** | Chinese herbology | TCM herbalism concept | ~200 |
| **Q165348** | acupuncture | Acupuncture therapy | ~100 |
| **Q16068429** | herbal formula | General herbal formulas | ~1,500 |
| **Q2993899** | kampo medicine | Japanese TCM | ~500 |

### 1.4 Ayurvedic Medicine

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q132325** | Ayurveda | Ayurvedic system | ~500 |
| **Q2629892** | Ayurvedic medicine | Ayurvedic preparations | ~2,000 |
| **Q854936** | Unani medicine | Unani system | ~200 |
| **Q115089** | homeopathic remedy | Homeopathic preparations | ~3,000 |
| **Q54896073** | traditional medicine | General traditional medicine | ~5,000 |
| **Q12094395** | herbal medicine | General herbal medicine | ~10,000 |

### 1.5 Kampo Medicine (Japanese)

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q2993899** | kampo | Kampo medicine system | ~500 |
| **Q51128287** | traditional Chinese medicine preparation | Includes Kampo | ~5,000 |
| **Q27136935** | Japanese Pharmacopoeia drug | JP listed drugs | ~2,500 |

### 1.6 Western Herbs and Botanicals

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q2225** | herb | Herbaceous plants | ~50,000 |
| **Q756** | plant | Plant kingdom (too broad) | ~400,000+ |
| **Q188617** | medicinal plant | Plants used medicinally | ~25,000 |
| **Q2596997** | phytomedicine | Plant-based medicine | ~3,000 |
| **Q114657** | phytotherapy | Herbal therapy concept | ~200 |
| **Q43656** | essential oil | Plant essential oils | ~5,000 |
| **Q33685** | spice | Culinary/medicinal spices | ~3,000 |
| **Q2095** | food | Food items (broad) | ~100,000+ |

### 1.7 Supplements, Vitamins, and Minerals

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q324546** | dietary supplement | Supplements | ~15,000 |
| **Q34956** | vitamin | Vitamins | ~100 |
| **Q12674** | mineral (nutrient) | Nutritional minerals | ~50 |
| **Q407071** | micronutrient | Micronutrients | ~200 |
| **Q44432** | nutrient | Nutrients general | ~2,000 |
| **Q131436** | amino acid | Amino acids | ~500 |
| **Q47512** | protein supplement | Protein supplements | ~300 |
| **Q422212** | probiotic | Probiotic products | ~500 |
| **Q1571780** | nutraceutical | Nutraceuticals | ~2,000 |
| **Q898273** | prebiotic | Prebiotic substances | ~100 |

### 1.8 Genes

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q7187** | gene | Gene (primary class) | ~3,000,000 |
| **Q20747295** | protein-coding gene | Protein-coding genes | ~500,000 |
| **Q277338** | pseudogene | Pseudogenes | ~50,000 |
| **Q427087** | non-coding RNA | ncRNA genes | ~100,000 |
| **Q7949102** | microRNA gene | miRNA genes | ~5,000 |
| **Q284578** | transfer RNA | tRNA genes | ~10,000 |
| **Q11053** | ribosomal RNA | rRNA genes | ~5,000 |
| **Q22269212** | protein-coding gene in humans | Human protein-coding | ~20,000 |
| **Q106345675** | human gene | Human genes all types | ~60,000 |

### 1.9 Proteins and Enzymes

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q8054** | protein | Protein (primary class) | ~5,000,000 |
| **Q8047** | enzyme | Enzymes | ~300,000 |
| **Q417841** | receptor | Receptor proteins | ~50,000 |
| **Q422073** | transporter protein | Transporters | ~30,000 |
| **Q68619** | transcription factor | TFs | ~50,000 |
| **Q185583** | kinase | Kinase enzymes | ~20,000 |
| **Q131029** | protease | Proteolytic enzymes | ~15,000 |
| **Q182897** | ion channel | Ion channels | ~10,000 |
| **Q898362** | cytochrome P450 | CYP enzymes | ~2,000 |
| **Q420927** | antibody | Antibody proteins | ~10,000 |

### 1.10 Pathways and Biological Processes

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q4915012** | biological pathway | Pathway (primary) | ~50,000 |
| **Q2996394** | biological process | GO:BP annotations | ~30,000 |
| **Q14860489** | metabolic pathway | Metabolic pathways | ~5,000 |
| **Q188907** | signal transduction | Signaling pathways | ~3,000 |
| **Q842908** | gene regulatory network | Regulatory networks | ~1,000 |
| **Q178593** | cell signaling | Cell signaling | ~2,000 |
| **Q7182** | gene expression | Expression processes | ~1,000 |
| **Q66055** | biochemical cascade | Cascades | ~500 |
| **Q44948** | apoptosis | Apoptosis pathway | ~200 |
| **Q128406** | inflammation | Inflammatory processes | ~500 |

### 1.11 Diseases and Conditions

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q12136** | disease | Disease (primary class) | ~200,000 |
| **Q929833** | rare disease | Rare diseases | ~10,000 |
| **Q18553442** | genetic disorder | Genetic diseases | ~8,000 |
| **Q169872** | cancer | Cancers | ~5,000 |
| **Q181257** | infectious disease | Infections | ~10,000 |
| **Q11651** | cardiovascular disease | CV diseases | ~2,000 |
| **Q12152** | autoimmune disease | Autoimmune conditions | ~1,500 |
| **Q60151** | neurological disorder | Neuro diseases | ~3,000 |
| **Q3286409** | mental disorder | Mental conditions | ~2,000 |
| **Q12206** | diabetes mellitus | Diabetes types | ~100 |
| **Q112193867** | disease phenotype | Phenotypes | ~50,000 |

### 1.12 Molecular Targets and Mechanisms

| Q-ID | Name | Description | Est. Items |
|------|------|-------------|------------|
| **Q417841** | receptor | Receptors (drug targets) | ~50,000 |
| **Q74224** | nuclear receptor | Nuclear receptors | ~500 |
| **Q288864** | G protein-coupled receptor | GPCRs | ~1,000 |
| **Q185583** | kinase | Kinase targets | ~20,000 |
| **Q68561** | epigenetic modification | Epigenetic marks | ~500 |
| **Q424902** | chromatin remodeling | Chromatin processes | ~200 |

---

## 2. Complete Property (P-code) Reference

### 2.1 Chemical Structure Identifiers

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P233** | SMILES | CC(=O)OC1=CC=CC=C1C(=O)O | Canonical SMILES |
| **P234** | InChI | InChI=1S/C9H8O4/c... | IUPAC InChI |
| **P235** | InChIKey | BSYNRYMUTXBXSQ-UHFFFAOYSA-N | InChI hash |
| **P231** | CAS Registry Number | 50-78-2 | CAS |
| **P2017** | isomeric SMILES | Stereochemistry | Isomeric |
| **P274** | chemical formula | C9H8O4 | Molecular formula |
| **P2054** | density | 1.4 g/cm3 | Physical property |
| **P2101** | melting point | 136 C | Physical property |
| **P2102** | boiling point | 140 C | Physical property |
| **P2128** | solubility | 3 g/L | Solubility |
| **P2064** | KOW | 1.19 | Partition coefficient |

### 2.2 Chemical Database Identifiers

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P662** | PubChem CID | 2244 | PubChem |
| **P661** | PubChem SID | 174073 | PubChem Substance |
| **P683** | ChEBI ID | 15365 | ChEBI |
| **P592** | ChEMBL ID | CHEMBL25 | ChEMBL |
| **P715** | DrugBank ID | DB00945 | DrugBank |
| **P595** | Guide to Pharmacology ID | 4139 | IUPHAR |
| **P2892** | UMLS CUI | C0004057 | UMLS |
| **P486** | MeSH descriptor ID | D001241 | MeSH |
| **P5270** | MeSH tree code | D03.633.100.064 | MeSH tree |
| **P652** | UNII | R16CO5Y76E | FDA UNII |
| **P2275** | WHO INN | aspirin | WHO name |
| **P3345** | RxNorm CUI | 1191 | RxNorm |
| **P3395** | ChemSpider ID | 2157 | ChemSpider |
| **P2566** | ECHA InfoCard ID | 100.000.059 | ECHA |
| **P2840** | HSDB ID | 168 | HSDB |
| **P3117** | DSSTox substance ID | DTXSID... | EPA DSSTox |
| **P2892** | UMLS CUI | C0004057 | UMLS |
| **P3636** | PDB ligand ID | ASP | PDB |
| **P8189** | ChEMBL target ID | CHEMBL203 | ChEMBL target |

### 2.3 Drug Classification Identifiers

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P267** | ATC code | N02BA01 | WHO ATC |
| **P3489** | pregnancy category | D | FDA/AU category |
| **P6680** | legal status (medicine) | OTC | Regulation |
| **P2868** | subject has role | Q427087 (enzyme inhibitor) | Pharmacological role |
| **P636** | route of administration | oral | Admin route |
| **P3780** | active ingredient in | Q... | Formulation link |
| **P3781** | has active ingredient | Q... | Ingredient link |
| **P2844** | absorbed by | small intestine | Absorption |
| **P3354** | positive therapeutic predictor | CYP2D6 | PGx predictor |
| **P3355** | negative therapeutic predictor | CYP2C19 | PGx predictor |
| **P769** | significant drug interaction | Q407548 | Drug interaction |

### 2.4 Biological Identifiers - Genes

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P351** | Entrez Gene ID | 7157 | NCBI Gene |
| **P353** | HGNC gene symbol | TP53 | HGNC |
| **P354** | HGNC ID | HGNC:11998 | HGNC |
| **P594** | Ensembl gene ID | ENSG00000141510 | Ensembl |
| **P593** | HomoloGene ID | 460 | HomoloGene |
| **P5806** | NCBI Taxonomy ID | 9606 | NCBI Taxonomy |
| **P2249** | RefSeq genome ID | NC_000017.11 | RefSeq |
| **P639** | RefSeq RNA ID | NM_000546.6 | RefSeq |
| **P637** | RefSeq protein ID | NP_000537.3 | RefSeq |
| **P2393** | NCBI Locus tag | TP53 | NCBI |
| **P704** | Ensembl transcript ID | ENST00000269305 | Ensembl |
| **P705** | Ensembl protein ID | ENSP00000269305 | Ensembl |
| **P684** | ortholog | Q... | Ortholog link |
| **P688** | encodes | Q... (protein) | Gene->Protein |
| **P702** | encoded by | Q... (gene) | Protein->Gene |
| **P703** | found in taxon | Q15978631 | Species |
| **P1057** | chromosome | 17 | Chromosome |
| **P644** | genomic start | 7661779 | Start position |
| **P645** | genomic end | 7687538 | End position |
| **P2548** | strand orientation | reverse strand | Strand |

### 2.5 Biological Identifiers - Proteins

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P352** | UniProt protein ID | P04637 | UniProt |
| **P705** | Ensembl protein ID | ENSP00000269305 | Ensembl |
| **P637** | RefSeq protein ID | NP_000537.3 | RefSeq |
| **P591** | EC enzyme number | 2.7.11.1 | EC |
| **P638** | PDB structure ID | 1TUP | PDB |
| **P3219** | Protein Atlas ID | ENSG00000141510 | HPA |
| **P680** | molecular function | GO:0003677 | GO MF |
| **P681** | cell component | GO:0005634 | GO CC |
| **P682** | biological process | GO:0006915 | GO BP |
| **P129** | physically interacts with | Q... | PPI |
| **P128** | regulates | Q... | Regulation |
| **P1910** | decreased expression in | Q... | Down in disease |
| **P1909** | increased expression in | Q... | Up in disease |

### 2.6 Pathway and Process Identifiers

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P3937** | Reactome pathway ID | R-HSA-109581 | Reactome |
| **P2888** | KEGG pathway ID | hsa04210 | KEGG |
| **P2410** | WikiPathways ID | WP254 | WikiPathways |
| **P680** | molecular function | GO:0003677 | GO |
| **P681** | cell component | GO:0005634 | GO |
| **P682** | biological process | GO:0006915 | GO |
| **P686** | Gene Ontology ID | GO:0006915 | GO |
| **P527** | has part | Q... | Pathway component |
| **P361** | part of | Q... | Part of pathway |

### 2.7 Medical and Clinical Properties

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P2175** | medical condition treated | Q41861 | Indication |
| **P1050** | medical condition | Q12136 | Disease link |
| **P780** | symptoms and signs | Q485831 | Symptoms |
| **P1995** | health specialty | Q7867 | Medical specialty |
| **P2176** | drug or therapy used for treatment | Q188724 | Treatment |
| **P1582** | natural product of taxon | Q7469 | Natural source |
| **P5572** | has cause | Q12136 | Disease cause |
| **P780** | symptoms and signs | Q485831 | Symptoms |
| **P1554** | risk factor | Q12174 | Risk factor |
| **P828** | has cause | Q... | Causation |
| **P1060** | disease transmission process | Q... | Transmission |

### 2.8 Pharmacological Relationships

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P129** | physically interacts with | Q412415 (COX-2) | Target binding |
| **P769** | significant drug interaction | Q407548 (warfarin) | DDI |
| **P3489** | mechanism of action | Q... | MOA |
| **P2868** | subject has role | Q427087 | Drug role |
| **P3354** | positive therapeutic predictor | Q... | PGx +ve |
| **P3355** | negative therapeutic predictor | Q... | PGx -ve |
| **P3356** | positive diagnostic predictor | Q... | Dx predictor |
| **P3433** | biological variant of | Q... | Variant |
| **P3364** | inhibitor of | Q... | Inhibition |
| **P3771** | activator of | Q... | Activation |
| **P4044** | therapeutic area | Q... | Therapy area |
| **P5167** | substrate of | Q... | Substrate |

### 2.9 Traditional Medicine Properties

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P2175** | medical condition treated | Q... | Indication |
| **P1582** | natural product of taxon | Q... | Plant source |
| **P927** | anatomical location | Q... | Body location |
| **P2176** | drug or therapy used | Q... | Therapy |
| **P1995** | health specialty | Q891104 | TCM specialty |
| **P279** | subclass of | Q891104 | TCM subclass |
| **P361** | part of | Q51128287 | Part of formula |
| **P527** | has part | Q... | Formula ingredient |
| **P6366** | Microsoft Academic ID | varies | Academic ID |

### 2.10 Publication and Evidence

| P-code | Name | Example Value | Database |
|--------|------|---------------|----------|
| **P698** | PubMed ID | 12505355 | PubMed |
| **P932** | PMC ID | PMC123456 | PMC |
| **P356** | DOI | 10.1038/nature12831 | DOI |
| **P248** | stated in | Q... | Source |
| **P813** | retrieved | 2024-01-15 | Date |
| **P854** | reference URL | https://... | URL |
| **P1433** | published in | Q... | Journal |
| **P50** | author | Q... | Author |
| **P577** | publication date | 2024-01-01 | Date |

### 2.11 Administrative and Classification

| P-code | Name | Description |
|--------|------|-------------|
| **P31** | instance of | Primary classification |
| **P279** | subclass of | Class hierarchy |
| **P361** | part of | Membership |
| **P527** | has part | Components |
| **P910** | topic's main category | Category |
| **P1269** | facet of | Aspect |
| **P460** | said to be the same as | Identity |
| **P1889** | different from | Disambiguation |
| **P2860** | cites work | Citation |

---

## 3. Dump Filtering Strategy

### 3.1 Overview

The Wikidata JSON dump is approximately 100GB compressed (1TB uncompressed). Efficient filtering requires:
1. Streaming processing (line-by-line)
2. Early termination based on P31 (instance of)
3. Property whitelist extraction

### 3.2 Filter by Instance-of (P31) Values

```python
#!/usr/bin/env python3
"""
Wikidata Dump Filter - Extract biomedical entities
"""

import bz2
import gzip
import json
import sys
from typing import Set, Dict, Iterator, Any
from pathlib import Path
from dataclasses import dataclass, field

# Master Q-ID sets for filtering
DRUG_QIDS = {
    'Q12140',      # medication
    'Q35456',      # pharmaceutical drug
    'Q8386',       # drug
    'Q422248',     # monoclonal antibody
    'Q17143904',   # small molecule drug
    'Q27136782',   # biological medicine
    'Q18214535',   # prodrug
    'Q79460',      # vaccine
}

CHEMICAL_QIDS = {
    'Q11173',      # chemical compound
    'Q79529',      # chemical substance
    'Q11344',      # alkaloid
    'Q134219',     # flavonoid
    'Q190875',     # terpenoid
    'Q407595',     # glycoside
    'Q193893',     # steroid
}

TCM_QIDS = {
    'Q891104',     # traditional Chinese medicine
    'Q51128287',   # TCM preparation
    'Q18025426',   # Chinese herbal medicine
    'Q54951330',   # TCM component
    'Q2993899',    # kampo medicine
}

AYURVEDA_QIDS = {
    'Q132325',     # Ayurveda
    'Q2629892',    # Ayurvedic medicine
    'Q854936',     # Unani medicine
    'Q115089',     # homeopathic remedy
}

HERB_QIDS = {
    'Q188617',     # medicinal plant
    'Q2596997',    # phytomedicine
    'Q43656',      # essential oil
    'Q12094395',   # herbal medicine
}

SUPPLEMENT_QIDS = {
    'Q324546',     # dietary supplement
    'Q34956',      # vitamin
    'Q12674',      # mineral nutrient
    'Q131436',     # amino acid
    'Q1571780',    # nutraceutical
}

GENE_QIDS = {
    'Q7187',       # gene
    'Q20747295',   # protein-coding gene
    'Q277338',     # pseudogene
    'Q427087',     # non-coding RNA
    'Q22269212',   # human protein-coding gene
}

PROTEIN_QIDS = {
    'Q8054',       # protein
    'Q8047',       # enzyme
    'Q417841',     # receptor
    'Q422073',     # transporter
    'Q68619',      # transcription factor
    'Q185583',     # kinase
    'Q898362',     # cytochrome P450
}

PATHWAY_QIDS = {
    'Q4915012',    # biological pathway
    'Q2996394',    # biological process
    'Q14860489',   # metabolic pathway
    'Q188907',     # signal transduction
}

DISEASE_QIDS = {
    'Q12136',      # disease
    'Q929833',     # rare disease
    'Q18553442',   # genetic disorder
    'Q169872',     # cancer
    'Q181257',     # infectious disease
}

# Combine all target Q-IDs
ALL_TARGET_QIDS = (
    DRUG_QIDS | CHEMICAL_QIDS | TCM_QIDS | AYURVEDA_QIDS |
    HERB_QIDS | SUPPLEMENT_QIDS | GENE_QIDS | PROTEIN_QIDS |
    PATHWAY_QIDS | DISEASE_QIDS
)

# Properties to extract
PROPERTIES_TO_EXTRACT = {
    # Chemical identifiers
    'P233', 'P234', 'P235', 'P231', 'P274', 'P662', 'P683',
    'P592', 'P715', 'P652', 'P267', 'P486', 'P2892', 'P3345',
    'P595', 'P2275', 'P3395',

    # Gene identifiers
    'P351', 'P353', 'P354', 'P594', 'P639', 'P637', 'P704', 'P705',

    # Protein identifiers
    'P352', 'P591', 'P638',

    # Pathway identifiers
    'P3937', 'P2888', 'P2410', 'P686', 'P680', 'P681', 'P682',

    # Medical properties
    'P2175', 'P769', 'P129', 'P2868', 'P3489', 'P636', 'P780',

    # Structural relationships
    'P31', 'P279', 'P361', 'P527', 'P688', 'P702', 'P703',

    # Traditional medicine
    'P1582', 'P927', 'P2176', 'P1995',

    # References
    'P698', 'P932', 'P356',
}


@dataclass
class FilterStats:
    """Statistics for filtering run."""
    total_processed: int = 0
    drugs_found: int = 0
    chemicals_found: int = 0
    tcm_found: int = 0
    genes_found: int = 0
    proteins_found: int = 0
    diseases_found: int = 0
    pathways_found: int = 0
    supplements_found: int = 0
    other_found: int = 0


def open_dump(filepath: str):
    """Open compressed or uncompressed dump file."""
    if filepath.endswith('.bz2'):
        return bz2.open(filepath, 'rt', encoding='utf-8')
    elif filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt', encoding='utf-8')
    else:
        return open(filepath, 'r', encoding='utf-8')


def get_instance_of_qids(entity: Dict) -> Set[str]:
    """Extract all Q-IDs from P31 (instance of) claims."""
    qids = set()
    claims = entity.get('claims', {})

    for claim in claims.get('P31', []):
        try:
            qid = claim['mainsnak']['datavalue']['value']['id']
            qids.add(qid)
        except (KeyError, TypeError):
            continue

    return qids


def get_subclass_of_qids(entity: Dict) -> Set[str]:
    """Extract all Q-IDs from P279 (subclass of) claims."""
    qids = set()
    claims = entity.get('claims', {})

    for claim in claims.get('P279', []):
        try:
            qid = claim['mainsnak']['datavalue']['value']['id']
            qids.add(qid)
        except (KeyError, TypeError):
            continue

    return qids


def classify_entity(instance_qids: Set[str], subclass_qids: Set[str]) -> str:
    """Classify entity into category based on Q-IDs."""
    all_qids = instance_qids | subclass_qids

    if all_qids & DRUG_QIDS:
        return 'drug'
    elif all_qids & GENE_QIDS:
        return 'gene'
    elif all_qids & PROTEIN_QIDS:
        return 'protein'
    elif all_qids & DISEASE_QIDS:
        return 'disease'
    elif all_qids & PATHWAY_QIDS:
        return 'pathway'
    elif all_qids & TCM_QIDS:
        return 'tcm'
    elif all_qids & AYURVEDA_QIDS:
        return 'ayurveda'
    elif all_qids & HERB_QIDS:
        return 'herb'
    elif all_qids & SUPPLEMENT_QIDS:
        return 'supplement'
    elif all_qids & CHEMICAL_QIDS:
        return 'chemical'
    else:
        return 'other'


def extract_properties(entity: Dict, properties: Set[str]) -> Dict:
    """Extract only specified properties from entity."""
    extracted = {
        'id': entity.get('id'),
        'type': entity.get('type'),
        'labels': {},
        'descriptions': {},
        'aliases': {},
        'claims': {},
        'sitelinks': {},
    }

    # Extract English labels/descriptions/aliases
    if 'labels' in entity and 'en' in entity['labels']:
        extracted['labels']['en'] = entity['labels']['en']

    if 'descriptions' in entity and 'en' in entity['descriptions']:
        extracted['descriptions']['en'] = entity['descriptions']['en']

    if 'aliases' in entity and 'en' in entity['aliases']:
        extracted['aliases']['en'] = entity['aliases']['en']

    # Extract only specified properties
    claims = entity.get('claims', {})
    for prop in properties:
        if prop in claims:
            extracted['claims'][prop] = claims[prop]

    # Extract Wikipedia sitelinks
    sitelinks = entity.get('sitelinks', {})
    if 'enwiki' in sitelinks:
        extracted['sitelinks']['enwiki'] = sitelinks['enwiki']

    return extracted


def stream_wikidata_dump(filepath: str) -> Iterator[Dict]:
    """Stream entities from Wikidata JSON dump."""
    with open_dump(filepath) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip JSON array brackets
            if line in ['[', ']', '']:
                continue

            # Remove trailing comma
            if line.endswith(','):
                line = line[:-1]

            try:
                entity = json.loads(line)
                yield entity
            except json.JSONDecodeError as e:
                if line_num % 1000000 == 0:
                    print(f"Warning: JSON error at line {line_num}: {e}",
                          file=sys.stderr)
                continue


def filter_wikidata_dump(
    input_path: str,
    output_path: str,
    target_qids: Set[str] = ALL_TARGET_QIDS,
    properties: Set[str] = PROPERTIES_TO_EXTRACT,
    progress_interval: int = 100000
) -> FilterStats:
    """
    Filter Wikidata dump for relevant entities.

    Args:
        input_path: Path to Wikidata JSON dump (bz2/gz/json)
        output_path: Path for filtered output (JSON lines)
        target_qids: Set of Q-IDs to filter for
        properties: Set of P-codes to extract
        progress_interval: Print progress every N entities

    Returns:
        FilterStats with counts
    """
    stats = FilterStats()

    with open(output_path, 'w', encoding='utf-8') as out_f:
        for entity in stream_wikidata_dump(input_path):
            stats.total_processed += 1

            # Progress reporting
            if stats.total_processed % progress_interval == 0:
                print(f"Processed: {stats.total_processed:,} | "
                      f"Found: drugs={stats.drugs_found}, genes={stats.genes_found}, "
                      f"proteins={stats.proteins_found}, diseases={stats.diseases_found}")

            # Skip non-items (properties, lexemes)
            if entity.get('type') != 'item':
                continue

            # Check P31 (instance of) and P279 (subclass of)
            instance_qids = get_instance_of_qids(entity)
            subclass_qids = get_subclass_of_qids(entity)
            all_qids = instance_qids | subclass_qids

            # Check if entity matches any target Q-ID
            if not (all_qids & target_qids):
                continue

            # Classify and count
            category = classify_entity(instance_qids, subclass_qids)

            if category == 'drug':
                stats.drugs_found += 1
            elif category == 'gene':
                stats.genes_found += 1
            elif category == 'protein':
                stats.proteins_found += 1
            elif category == 'disease':
                stats.diseases_found += 1
            elif category == 'pathway':
                stats.pathways_found += 1
            elif category == 'tcm':
                stats.tcm_found += 1
            elif category == 'chemical':
                stats.chemicals_found += 1
            elif category == 'supplement':
                stats.supplements_found += 1
            else:
                stats.other_found += 1

            # Extract relevant properties and save
            filtered_entity = extract_properties(entity, properties)
            filtered_entity['_category'] = category

            out_f.write(json.dumps(filtered_entity, ensure_ascii=False) + '\n')

    return stats


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Filter Wikidata dump for biomedical entities'
    )
    parser.add_argument('input', help='Input Wikidata dump path')
    parser.add_argument('output', help='Output JSON lines path')
    parser.add_argument('--progress', type=int, default=100000,
                        help='Progress interval')

    args = parser.parse_args()

    print(f"Filtering {args.input} -> {args.output}")
    print(f"Target Q-IDs: {len(ALL_TARGET_QIDS)}")
    print(f"Properties to extract: {len(PROPERTIES_TO_EXTRACT)}")

    stats = filter_wikidata_dump(
        args.input,
        args.output,
        progress_interval=args.progress
    )

    print("\n=== Final Statistics ===")
    print(f"Total processed: {stats.total_processed:,}")
    print(f"Drugs: {stats.drugs_found:,}")
    print(f"Genes: {stats.genes_found:,}")
    print(f"Proteins: {stats.proteins_found:,}")
    print(f"Diseases: {stats.diseases_found:,}")
    print(f"Pathways: {stats.pathways_found:,}")
    print(f"TCM: {stats.tcm_found:,}")
    print(f"Chemicals: {stats.chemicals_found:,}")
    print(f"Supplements: {stats.supplements_found:,}")
    print(f"Other: {stats.other_found:,}")
```

### 3.3 Output Formats

#### JSON Lines (Recommended for streaming)
```bash
# Each line is a complete JSON object
{"id":"Q188724","labels":{"en":{"value":"Ibuprofen"}},"claims":{...}}
{"id":"Q407431","labels":{"en":{"value":"Metformin"}},"claims":{...}}
```

#### TSV Export Script
```python
def export_to_tsv(jsonl_path: str, tsv_path: str):
    """Convert filtered JSON lines to TSV."""
    import csv

    with open(jsonl_path, 'r') as f_in, \
         open(tsv_path, 'w', newline='') as f_out:

        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow([
            'qid', 'label', 'category', 'drugbank_id', 'pubchem_cid',
            'chebi_id', 'chembl_id', 'uniprot_id', 'entrez_id', 'smiles'
        ])

        for line in f_in:
            entity = json.loads(line)

            # Extract values
            qid = entity['id']
            label = entity.get('labels', {}).get('en', {}).get('value', '')
            category = entity.get('_category', '')

            claims = entity.get('claims', {})

            def get_claim_value(prop):
                if prop not in claims:
                    return ''
                try:
                    return claims[prop][0]['mainsnak']['datavalue']['value']
                except (KeyError, IndexError, TypeError):
                    return ''

            def get_claim_id(prop):
                val = get_claim_value(prop)
                if isinstance(val, dict):
                    return val.get('id', '')
                return str(val) if val else ''

            writer.writerow([
                qid,
                label,
                category,
                get_claim_id('P715'),   # DrugBank
                get_claim_id('P662'),   # PubChem
                get_claim_id('P683'),   # ChEBI
                get_claim_id('P592'),   # ChEMBL
                get_claim_id('P352'),   # UniProt
                get_claim_id('P351'),   # Entrez
                get_claim_id('P233'),   # SMILES
            ])
```

#### SQLite Export Script
```python
def export_to_sqlite(jsonl_path: str, db_path: str):
    """Convert filtered JSON lines to SQLite database."""
    import sqlite3

    conn = sqlite3.connect(db_path)
    conn.execute('''
        CREATE TABLE IF NOT EXISTS entities (
            qid TEXT PRIMARY KEY,
            label TEXT,
            description TEXT,
            category TEXT,
            claims_json TEXT
        )
    ''')

    conn.execute('''
        CREATE TABLE IF NOT EXISTS identifiers (
            qid TEXT,
            property TEXT,
            value TEXT,
            FOREIGN KEY (qid) REFERENCES entities(qid)
        )
    ''')

    conn.execute('CREATE INDEX IF NOT EXISTS idx_identifiers_prop ON identifiers(property)')
    conn.execute('CREATE INDEX IF NOT EXISTS idx_identifiers_val ON identifiers(value)')

    with open(jsonl_path, 'r') as f:
        batch_entities = []
        batch_identifiers = []

        for line in f:
            entity = json.loads(line)

            qid = entity['id']
            label = entity.get('labels', {}).get('en', {}).get('value', '')
            desc = entity.get('descriptions', {}).get('en', {}).get('value', '')
            category = entity.get('_category', '')

            batch_entities.append((qid, label, desc, category, json.dumps(entity['claims'])))

            # Extract identifier properties
            for prop in ['P715', 'P662', 'P683', 'P592', 'P352', 'P351', 'P233',
                         'P234', 'P235', 'P231', 'P267', 'P353', 'P354', 'P594']:
                claims = entity.get('claims', {}).get(prop, [])
                for claim in claims:
                    try:
                        val = claim['mainsnak']['datavalue']['value']
                        if isinstance(val, dict):
                            val = val.get('id', str(val))
                        batch_identifiers.append((qid, prop, str(val)))
                    except (KeyError, TypeError):
                        continue

            if len(batch_entities) >= 10000:
                conn.executemany(
                    'INSERT OR REPLACE INTO entities VALUES (?,?,?,?,?)',
                    batch_entities
                )
                conn.executemany(
                    'INSERT INTO identifiers VALUES (?,?,?)',
                    batch_identifiers
                )
                conn.commit()
                batch_entities = []
                batch_identifiers = []

        # Final batch
        if batch_entities:
            conn.executemany(
                'INSERT OR REPLACE INTO entities VALUES (?,?,?,?,?)',
                batch_entities
            )
            conn.executemany(
                'INSERT INTO identifiers VALUES (?,?,?)',
                batch_identifiers
            )
            conn.commit()

    conn.close()
```

---

## 4. Master SPARQL Queries

### 4.1 All Drugs with All Identifiers

```sparql
# Query: Get all pharmaceutical drugs with complete identifier set
# Endpoint: https://query.wikidata.org/sparql
# Timeout: May need pagination for full dataset

SELECT DISTINCT
  ?drug ?drugLabel ?drugDescription
  ?cas ?drugbank ?pubchem ?chebi ?chembl ?unii ?atc ?mesh
  ?smiles ?inchi ?inchikey ?formula
  ?target ?targetLabel ?targetUniProt
  ?indication ?indicationLabel
WHERE {
  # Base: pharmaceutical drugs
  ?drug wdt:P31/wdt:P279* wd:Q12140 .

  # Chemical identifiers (all optional)
  OPTIONAL { ?drug wdt:P231 ?cas }
  OPTIONAL { ?drug wdt:P715 ?drugbank }
  OPTIONAL { ?drug wdt:P662 ?pubchem }
  OPTIONAL { ?drug wdt:P683 ?chebi }
  OPTIONAL { ?drug wdt:P592 ?chembl }
  OPTIONAL { ?drug wdt:P652 ?unii }
  OPTIONAL { ?drug wdt:P267 ?atc }
  OPTIONAL { ?drug wdt:P486 ?mesh }

  # Structure
  OPTIONAL { ?drug wdt:P233 ?smiles }
  OPTIONAL { ?drug wdt:P234 ?inchi }
  OPTIONAL { ?drug wdt:P235 ?inchikey }
  OPTIONAL { ?drug wdt:P274 ?formula }

  # Targets
  OPTIONAL {
    ?drug wdt:P129 ?target .
    ?target wdt:P31 wd:Q8054 .
    OPTIONAL { ?target wdt:P352 ?targetUniProt }
  }

  # Indications
  OPTIONAL { ?drug wdt:P2175 ?indication }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 50000
```

### 4.2 All Genes with All Identifiers

```sparql
# Query: Get all human genes with complete identifier set
# Note: Filter to human genes (Q15978631) for manageable results

SELECT DISTINCT
  ?gene ?geneLabel ?geneDescription
  ?entrez ?hgnc_symbol ?hgnc_id ?ensembl
  ?refseq_rna ?refseq_protein
  ?chromosome ?start ?end ?strand
  ?protein ?proteinLabel ?uniprot
  ?goMF ?goMFLabel
  ?goBP ?goBPLabel
WHERE {
  # Human protein-coding genes
  ?gene wdt:P31 wd:Q20747295 ;
        wdt:P703 wd:Q15978631 .  # Homo sapiens

  # Gene identifiers
  OPTIONAL { ?gene wdt:P351 ?entrez }
  OPTIONAL { ?gene wdt:P353 ?hgnc_symbol }
  OPTIONAL { ?gene wdt:P354 ?hgnc_id }
  OPTIONAL { ?gene wdt:P594 ?ensembl }
  OPTIONAL { ?gene wdt:P639 ?refseq_rna }
  OPTIONAL { ?gene wdt:P637 ?refseq_protein }

  # Genomic location
  OPTIONAL { ?gene wdt:P1057 ?chromosome }
  OPTIONAL { ?gene wdt:P644 ?start }
  OPTIONAL { ?gene wdt:P645 ?end }
  OPTIONAL { ?gene wdt:P2548 ?strand }

  # Encoded protein
  OPTIONAL {
    ?gene wdt:P688 ?protein .
    OPTIONAL { ?protein wdt:P352 ?uniprot }
  }

  # GO annotations (sample)
  OPTIONAL { ?gene wdt:P680 ?goMF }
  OPTIONAL { ?gene wdt:P682 ?goBP }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 30000
```

### 4.3 All Drug-Gene Interactions

```sparql
# Query: Drug-Target-Gene relationships
# Links drugs to their molecular targets and encoding genes

SELECT DISTINCT
  ?drug ?drugLabel ?drugbank
  ?target ?targetLabel ?targetUniProt
  ?gene ?geneLabel ?geneSymbol ?entrez
  ?interactionType ?interactionTypeLabel
WHERE {
  # Drugs with targets
  ?drug wdt:P31/wdt:P279* wd:Q12140 .

  # Drug identifiers
  OPTIONAL { ?drug wdt:P715 ?drugbank }

  # Drug-protein interaction
  ?drug wdt:P129 ?target .
  ?target wdt:P31 wd:Q8054 .

  # Protein identifiers
  OPTIONAL { ?target wdt:P352 ?targetUniProt }

  # Protein encoded by gene
  OPTIONAL {
    ?target wdt:P702 ?gene .
    ?gene wdt:P703 wd:Q15978631 .  # Human genes
    OPTIONAL { ?gene wdt:P353 ?geneSymbol }
    OPTIONAL { ?gene wdt:P351 ?entrez }
  }

  # Interaction type qualifier
  OPTIONAL {
    ?drug p:P129 ?interactionStatement .
    ?interactionStatement ps:P129 ?target .
    OPTIONAL { ?interactionStatement pq:P2868 ?interactionType }
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 100000
```

### 4.4 All Traditional Medicine Items

```sparql
# Query: TCM, Ayurveda, Kampo, and other traditional medicines

SELECT DISTINCT
  ?item ?itemLabel ?itemDescription
  ?system ?systemLabel
  ?source_plant ?source_plantLabel
  ?indication ?indicationLabel
  ?compound ?compoundLabel ?pubchem
WHERE {
  # Traditional medicine classes
  VALUES ?class {
    wd:Q891104      # Traditional Chinese Medicine
    wd:Q51128287    # TCM preparation
    wd:Q18025426    # Chinese herbal medicine
    wd:Q2993899     # Kampo
    wd:Q132325      # Ayurveda
    wd:Q2629892     # Ayurvedic medicine
    wd:Q12094395    # herbal medicine
    wd:Q2596997     # phytomedicine
    wd:Q188617      # medicinal plant
  }

  ?item wdt:P31/wdt:P279* ?class .

  # Traditional system
  OPTIONAL { ?item wdt:P1995 ?system }

  # Natural source
  OPTIONAL { ?item wdt:P1582 ?source_plant }

  # Indications
  OPTIONAL { ?item wdt:P2175 ?indication }

  # Active compounds
  OPTIONAL {
    ?item wdt:P527 ?compound .
    ?compound wdt:P31 wd:Q11173 .
    OPTIONAL { ?compound wdt:P662 ?pubchem }
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 50000
```

### 4.5 Federated Query - Wikidata + WikiPathways

```sparql
# Query: Link Wikidata drugs to WikiPathways pathways via ChEBI
# Note: Federated queries may timeout - use smaller batches

PREFIX wp: <http://vocabularies.wikipathways.org/wp#>
PREFIX dcterms: <http://purl.org/dc/terms/>

SELECT DISTINCT
  ?drug ?drugLabel ?chebiId
  ?pathwayId ?pathwayTitle
  ?gene ?geneLabel
WHERE {
  # Get drugs with ChEBI IDs from Wikidata
  ?drug wdt:P31 wd:Q12140 ;
        wdt:P683 ?chebiId .

  BIND(IRI(CONCAT("http://identifiers.org/chebi/CHEBI:", ?chebiId)) AS ?chebiUri)

  # Federated query to WikiPathways
  SERVICE <https://sparql.wikipathways.org/sparql> {
    ?metabolite a wp:Metabolite ;
                wp:bdbChEBI ?chebiUri ;
                dcterms:isPartOf ?pathway .
    ?pathway a wp:Pathway ;
             dcterms:identifier ?pathwayId ;
             dcterms:title ?pathwayTitle ;
             wp:organismName "Homo sapiens" .

    # Get genes in same pathway
    OPTIONAL {
      ?geneNode a wp:GeneProduct ;
                dcterms:isPartOf ?pathway ;
                rdfs:label ?gene .
    }
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 10000
```

### 4.6 Federated Query - Wikidata + Reactome

```sparql
# Query: Link Wikidata proteins to Reactome pathways

SELECT DISTINCT
  ?protein ?proteinLabel ?uniprot
  ?reactomeId ?pathwayName
WHERE {
  # Get proteins with UniProt IDs
  ?protein wdt:P31 wd:Q8054 ;
           wdt:P352 ?uniprot ;
           wdt:P703 wd:Q15978631 .  # Human

  # Get Reactome pathway
  OPTIONAL { ?protein wdt:P3937 ?reactomeId }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }

  # Note: For full Reactome data, query Reactome API separately
  BIND(
    IF(BOUND(?reactomeId),
       IRI(CONCAT("https://reactome.org/content/detail/", ?reactomeId)),
       ?nothing)
    AS ?reactomeLink
  )
}
LIMIT 20000
```

### 4.7 Pharmacogenomics Relationships

```sparql
# Query: Drug-Gene pharmacogenomic relationships

SELECT DISTINCT
  ?drug ?drugLabel ?drugbank
  ?gene ?geneLabel ?geneSymbol
  ?predictorType
  ?variant ?variantLabel
WHERE {
  # Drugs with PGx relationships
  ?drug wdt:P31 wd:Q12140 .
  OPTIONAL { ?drug wdt:P715 ?drugbank }

  # Positive predictor genes
  {
    ?drug wdt:P3354 ?gene .
    BIND("positive_predictor" AS ?predictorType)
  }
  UNION
  # Negative predictor genes
  {
    ?drug wdt:P3355 ?gene .
    BIND("negative_predictor" AS ?predictorType)
  }

  ?gene wdt:P703 wd:Q15978631 .  # Human
  OPTIONAL { ?gene wdt:P353 ?geneSymbol }

  # Associated variants
  OPTIONAL {
    ?drug wdt:P3433 ?variant .
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 10000
```

### 4.8 Disease-Gene-Drug Network

```sparql
# Query: Complete disease-gene-drug network

SELECT DISTINCT
  ?disease ?diseaseLabel ?mesh
  ?gene ?geneLabel ?geneSymbol
  ?drug ?drugLabel ?drugbank
  ?relationshipType
WHERE {
  # Diseases with gene associations
  ?disease wdt:P31/wdt:P279* wd:Q12136 .
  OPTIONAL { ?disease wdt:P486 ?mesh }

  # Disease-gene relationships
  {
    ?disease wdt:P2293 ?gene .  # genetic association
    BIND("genetic_association" AS ?relationshipType)
  }
  UNION
  {
    ?gene wdt:P1909 ?disease .  # increased expression
    BIND("increased_expression" AS ?relationshipType)
  }
  UNION
  {
    ?gene wdt:P1910 ?disease .  # decreased expression
    BIND("decreased_expression" AS ?relationshipType)
  }

  ?gene wdt:P703 wd:Q15978631 .  # Human
  OPTIONAL { ?gene wdt:P353 ?geneSymbol }

  # Drugs treating the disease
  OPTIONAL {
    ?drug wdt:P2175 ?disease .
    OPTIONAL { ?drug wdt:P715 ?drugbank }
  }

  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 50000
```

---

## 5. Data Completeness Matrix

### 5.1 Drug Data Coverage

| Data Element | Wikidata Est. | DrugBank | ChEMBL | Coverage vs DrugBank |
|--------------|---------------|----------|--------|----------------------|
| Total Drugs | ~45,000 | ~15,000 | ~2.4M | 300% (broader definition) |
| DrugBank ID (P715) | ~12,000 | 15,000 | - | 80% |
| PubChem CID (P662) | ~40,000 | ~13,000 | - | 300% |
| ChEMBL ID (P592) | ~30,000 | - | 2.4M | 1.3% |
| ATC Code (P267) | ~8,000 | ~5,500 | - | 145% |
| SMILES (P233) | ~35,000 | ~15,000 | ~2.2M | 230% vs DrugBank |
| CAS Number (P231) | ~25,000 | ~14,000 | - | 180% |
| Drug Targets (P129) | ~15,000 | ~19,000 | ~15M | 79% |
| Indications (P2175) | ~20,000 | ~6,000 | - | 330% |
| Drug-Drug Interactions | ~5,000 | ~300K | - | 1.7% |

### 5.2 Gene Data Coverage

| Data Element | Wikidata Est. | NCBI Gene | Ensembl | Coverage vs NCBI |
|--------------|---------------|-----------|---------|------------------|
| Human Genes | ~60,000 | ~60,000 | ~67,000 | 100% |
| Entrez ID (P351) | ~55,000 | 60,000 | - | 92% |
| HGNC Symbol (P353) | ~42,000 | ~42,000 | - | 100% |
| Ensembl ID (P594) | ~50,000 | - | 67,000 | 75% |
| GO Annotations | ~100,000 | ~500,000 | - | 20% |
| Gene-Disease | ~30,000 | ~20,000 | - | 150% |

### 5.3 Protein Data Coverage

| Data Element | Wikidata Est. | UniProt | PDB | Coverage vs UniProt |
|--------------|---------------|---------|-----|---------------------|
| Total Proteins | ~500,000 | ~250M | 220K | 0.2% (human-focused) |
| Human Proteins | ~20,000 | ~80,000 | - | 25% |
| UniProt ID (P352) | ~200,000 | 250M | - | 0.08% |
| PDB ID (P638) | ~15,000 | - | 220K | 7% |
| EC Number (P591) | ~8,000 | ~8,000 | - | 100% |

### 5.4 Traditional Medicine Coverage

| Data Element | Wikidata Est. | BATMAN-TCM | TCMBank | IMPPAT | Coverage |
|--------------|---------------|------------|---------|--------|----------|
| TCM Herbs | ~3,000 | 8,404 | 9,192 | - | 33-36% |
| TCM Compounds | ~5,000 | 39,171 | 61,966 | - | 8-13% |
| Ayurvedic Herbs | ~1,000 | - | - | 4,010 | 25% |
| Kampo Formulas | ~500 | - | - | - | Variable |

### 5.5 Pathway Data Coverage

| Data Element | Wikidata Est. | Reactome | WikiPathways | KEGG | Coverage |
|--------------|---------------|----------|--------------|------|----------|
| Pathways | ~50,000 | ~2,600 | ~3,000 | ~550 | 800-9000% |
| Reactome ID (P3937) | ~5,000 | 2,600 | - | - | 192% |
| WikiPathways ID (P2410) | ~2,000 | - | 3,000 | - | 67% |
| GO Process (P682) | ~30,000 | - | - | - | Variable |

### 5.6 Summary Statistics

| Category | Wikidata Unique Items | Primary Specialized DB | Wikidata Completeness |
|----------|----------------------|------------------------|----------------------|
| Drugs (approved) | ~45,000 | DrugBank: 15,000 | High (300%) |
| Chemical compounds | ~2.5M | PubChem: 115M | Low (2%) |
| Human genes | ~60,000 | NCBI Gene: 60,000 | High (100%) |
| Human proteins | ~20,000 | UniProt: 80,000 | Medium (25%) |
| TCM/Traditional | ~10,000 | BATMAN-TCM: 54,832 | Medium (18%) |
| Pathways | ~50,000 | Reactome: 2,600 | High (but less curated) |
| Diseases | ~200,000 | OMIM: 16,000 | High (1250%) |

### 5.7 Recommendations

**Use Wikidata For:**
- Cross-database identifier mapping (excellent)
- Drug-disease relationships (good coverage)
- Gene-disease relationships (good coverage)
- General entity discovery
- Traditional medicine starting point

**Supplement with Specialized DBs For:**
- Drug structures (ChEMBL/PubChem)
- Drug targets (ChEMBL/BindingDB)
- Protein details (UniProt)
- Pathways (Reactome/WikiPathways native)
- TCM compounds (BATMAN-TCM/TCMBank)
- Drug-drug interactions (DrugBank)

---

## 6. Recommended Extraction Pipeline

### 6.1 Pipeline Architecture

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     WIKIDATA EXTRACTION PIPELINE                             │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                   │
│  │   Wikidata   │───>│   Filter     │───>│   Extract    │                   │
│  │  JSON Dump   │    │   by Q-ID    │    │  Properties  │                   │
│  │  (~100 GB)   │    │              │    │              │                   │
│  └──────────────┘    └──────────────┘    └──────────────┘                   │
│         │                   │                   │                            │
│         │                   │                   │                            │
│         ▼                   ▼                   ▼                            │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐                   │
│  │   Download   │    │  Categories  │    │   Filtered   │                   │
│  │  with aria2c │    │  - drugs     │    │  JSON Lines  │                   │
│  │  (16 conn)   │    │  - genes     │    │  (~5 GB)     │                   │
│  └──────────────┘    │  - proteins  │    └──────────────┘                   │
│                      │  - diseases  │           │                            │
│                      │  - TCM       │           │                            │
│                      │  - pathways  │           ▼                            │
│                      └──────────────┘    ┌──────────────┐                   │
│                                          │   Normalize   │                   │
│                                          │   & Clean     │                   │
│                                          │              │                    │
│                                          └──────────────┘                   │
│                                                 │                            │
│                            ┌────────────────────┼────────────────────┐       │
│                            │                    │                    │       │
│                            ▼                    ▼                    ▼       │
│                     ┌──────────────┐    ┌──────────────┐    ┌──────────────┐│
│                     │   SQLite     │    │   TSV        │    │   ID Maps    ││
│                     │   Database   │    │   Exports    │    │   (JSON)     ││
│                     └──────────────┘    └──────────────┘    └──────────────┘│
│                            │                    │                    │       │
│                            └────────────────────┼────────────────────┘       │
│                                                 │                            │
│                                                 ▼                            │
│                                          ┌──────────────┐                   │
│                                          │   Merge with  │                   │
│                                          │  Specialized  │                   │
│                                          │     DBs       │                   │
│                                          └──────────────┘                   │
│                                                 │                            │
│                            ┌────────────────────┼────────────────────┐       │
│                            │                    │                    │       │
│                            ▼                    ▼                    ▼       │
│                     ┌──────────────┐    ┌──────────────┐    ┌──────────────┐│
│                     │   DrugBank   │    │   UniProt    │    │   Reactome   ││
│                     │   ChEMBL     │    │   NCBI Gene  │    │   WikiPW     ││
│                     │  BATMAN-TCM  │    │   Ensembl    │    │   KEGG       ││
│                     └──────────────┘    └──────────────┘    └──────────────┘│
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 6.2 Step-by-Step Commands

```bash
#!/bin/bash
# Complete Wikidata extraction pipeline

# 1. Download dump with multi-connection
echo "Step 1: Downloading Wikidata dump..."
aria2c -x 16 -s 16 \
    https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2 \
    -o wikidata-dump.json.bz2

# 2. Verify checksum
echo "Step 2: Verifying checksum..."
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2.md5
md5sum -c latest-all.json.bz2.md5

# 3. Filter for biomedical entities
echo "Step 3: Filtering for biomedical entities..."
python3 wikidata_filter.py \
    wikidata-dump.json.bz2 \
    wikidata-biomedical.jsonl \
    --progress 500000

# 4. Export to SQLite
echo "Step 4: Creating SQLite database..."
python3 wikidata_to_sqlite.py \
    wikidata-biomedical.jsonl \
    wikidata-biomedical.db

# 5. Create ID mapping files
echo "Step 5: Creating ID mapping files..."
python3 create_id_maps.py \
    wikidata-biomedical.jsonl \
    --output-dir ./id_maps/

# 6. Create TSV exports
echo "Step 6: Creating TSV exports..."
python3 export_tsv.py \
    wikidata-biomedical.jsonl \
    --output-dir ./tsv_exports/

echo "Pipeline complete!"
```

### 6.3 Resource Requirements

| Stage | CPU | RAM | Disk | Time |
|-------|-----|-----|------|------|
| Download | 1 core | 1 GB | 100 GB | 2-4 hours |
| Filter | 1-4 cores | 4-8 GB | 200 GB | 6-12 hours |
| SQLite Export | 1 core | 2-4 GB | 50 GB | 1-2 hours |
| TSV Export | 1 core | 2 GB | 10 GB | 30 min |
| **Total** | **4 cores** | **8 GB** | **400 GB** | **10-20 hours** |

---

## 7. Complete Wikidata Extractor Code

### 7.1 Main Extractor Script

```python
#!/usr/bin/env python3
"""
Complete Wikidata Extractor for Genetics/Health Knowledge Base

Features:
- Streams Wikidata JSON dump (handles 100GB+ files)
- Filters by Q-ID (instance-of and subclass-of)
- Extracts relevant properties only
- Outputs to multiple formats (JSON lines, SQLite, TSV)
- Includes progress tracking and resumption support

Usage:
    python wikidata_complete_extractor.py \
        --input latest-all.json.bz2 \
        --output-dir ./extracted/ \
        --formats jsonl sqlite tsv
"""

import argparse
import bz2
import gzip
import json
import sqlite3
import csv
import sys
import hashlib
from pathlib import Path
from typing import Set, Dict, Iterator, Any, List, Optional
from dataclasses import dataclass, field, asdict
from datetime import datetime
from collections import defaultdict
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class ExtractionConfig:
    """Configuration for Wikidata extraction."""

    # Q-IDs to filter for (instance-of or subclass-of)
    drug_qids: Set[str] = field(default_factory=lambda: {
        'Q12140', 'Q35456', 'Q8386', 'Q422248', 'Q17143904',
        'Q27136782', 'Q18214535', 'Q79460', 'Q58624061'
    })

    chemical_qids: Set[str] = field(default_factory=lambda: {
        'Q11173', 'Q79529', 'Q11344', 'Q134219', 'Q190875',
        'Q407595', 'Q2725376', 'Q193893', 'Q3314483'
    })

    tcm_qids: Set[str] = field(default_factory=lambda: {
        'Q891104', 'Q51128287', 'Q18025426', 'Q54951330',
        'Q899494', 'Q2993899', 'Q16068429'
    })

    ayurveda_qids: Set[str] = field(default_factory=lambda: {
        'Q132325', 'Q2629892', 'Q854936', 'Q115089', 'Q54896073'
    })

    herb_qids: Set[str] = field(default_factory=lambda: {
        'Q188617', 'Q2596997', 'Q43656', 'Q12094395', 'Q2225'
    })

    supplement_qids: Set[str] = field(default_factory=lambda: {
        'Q324546', 'Q34956', 'Q12674', 'Q407071', 'Q44432',
        'Q131436', 'Q47512', 'Q422212', 'Q1571780', 'Q898273'
    })

    gene_qids: Set[str] = field(default_factory=lambda: {
        'Q7187', 'Q20747295', 'Q277338', 'Q427087', 'Q7949102',
        'Q284578', 'Q11053', 'Q22269212', 'Q106345675'
    })

    protein_qids: Set[str] = field(default_factory=lambda: {
        'Q8054', 'Q8047', 'Q417841', 'Q422073', 'Q68619',
        'Q185583', 'Q131029', 'Q182897', 'Q898362', 'Q420927'
    })

    pathway_qids: Set[str] = field(default_factory=lambda: {
        'Q4915012', 'Q2996394', 'Q14860489', 'Q188907',
        'Q842908', 'Q178593', 'Q7182', 'Q66055'
    })

    disease_qids: Set[str] = field(default_factory=lambda: {
        'Q12136', 'Q929833', 'Q18553442', 'Q169872', 'Q181257',
        'Q11651', 'Q12152', 'Q60151', 'Q3286409', 'Q12206'
    })

    # Properties to extract
    properties: Set[str] = field(default_factory=lambda: {
        # Chemical identifiers
        'P233', 'P234', 'P235', 'P231', 'P274', 'P2017',
        'P662', 'P661', 'P683', 'P592', 'P715', 'P595',
        'P2892', 'P486', 'P5270', 'P652', 'P2275', 'P3345',
        'P3395', 'P2566', 'P2840', 'P3117', 'P3636', 'P8189',

        # Drug classification
        'P267', 'P3489', 'P6680', 'P2868', 'P636', 'P3780',
        'P3781', 'P2844', 'P3354', 'P3355', 'P769',

        # Gene identifiers
        'P351', 'P353', 'P354', 'P594', 'P593', 'P5806',
        'P2249', 'P639', 'P637', 'P2393', 'P704', 'P705',
        'P684', 'P688', 'P702', 'P703', 'P1057', 'P644',
        'P645', 'P2548',

        # Protein identifiers
        'P352', 'P591', 'P638', 'P3219',

        # GO and pathway
        'P680', 'P681', 'P682', 'P686', 'P3937', 'P2888',
        'P2410',

        # Medical
        'P2175', 'P1050', 'P780', 'P1995', 'P2176', 'P1582',
        'P5572', 'P1554', 'P828', 'P1060',

        # Pharmacology
        'P129', 'P128', 'P1910', 'P1909', 'P3364', 'P3771',
        'P4044', 'P5167', 'P3433',

        # Structure
        'P31', 'P279', 'P361', 'P527',

        # References
        'P698', 'P932', 'P356', 'P248', 'P813', 'P854',
    })

    @property
    def all_target_qids(self) -> Set[str]:
        """Get all Q-IDs to filter for."""
        return (
            self.drug_qids | self.chemical_qids | self.tcm_qids |
            self.ayurveda_qids | self.herb_qids | self.supplement_qids |
            self.gene_qids | self.protein_qids | self.pathway_qids |
            self.disease_qids
        )


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ExtractionStats:
    """Statistics for extraction run."""
    total_processed: int = 0
    total_matched: int = 0
    by_category: Dict[str, int] = field(default_factory=lambda: defaultdict(int))
    start_time: datetime = field(default_factory=datetime.now)

    def elapsed_time(self) -> float:
        return (datetime.now() - self.start_time).total_seconds()

    def processing_rate(self) -> float:
        elapsed = self.elapsed_time()
        return self.total_processed / elapsed if elapsed > 0 else 0


@dataclass
class ExtractedEntity:
    """Extracted and normalized entity."""
    qid: str
    label: str
    description: str
    aliases: List[str]
    category: str
    instance_of: List[str]
    subclass_of: List[str]
    identifiers: Dict[str, List[str]]
    relationships: Dict[str, List[Dict]]
    wikipedia_url: Optional[str]
    raw_claims: Dict


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def open_dump(filepath: str):
    """Open compressed or uncompressed dump file."""
    if filepath.endswith('.bz2'):
        return bz2.open(filepath, 'rt', encoding='utf-8')
    elif filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt', encoding='utf-8')
    else:
        return open(filepath, 'r', encoding='utf-8')


def stream_entities(filepath: str) -> Iterator[Dict]:
    """Stream entities from Wikidata JSON dump."""
    with open_dump(filepath) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Skip JSON array brackets
            if line in ['[', ']', '']:
                continue

            # Remove trailing comma
            if line.endswith(','):
                line = line[:-1]

            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue


def get_claim_values(claims: Dict, prop: str) -> List[Any]:
    """Extract all values for a property from claims."""
    values = []
    for claim in claims.get(prop, []):
        try:
            datavalue = claim.get('mainsnak', {}).get('datavalue', {})
            value = datavalue.get('value')

            if value is not None:
                # Handle different value types
                if isinstance(value, dict):
                    if 'id' in value:
                        values.append(value['id'])
                    elif 'amount' in value:
                        values.append(value['amount'])
                    elif 'text' in value:
                        values.append(value['text'])
                    elif 'time' in value:
                        values.append(value['time'])
                    else:
                        values.append(str(value))
                else:
                    values.append(value)
        except (KeyError, TypeError):
            continue

    return values


def get_claim_with_qualifiers(claims: Dict, prop: str) -> List[Dict]:
    """Extract claims with qualifiers for relationship properties."""
    results = []
    for claim in claims.get(prop, []):
        try:
            datavalue = claim.get('mainsnak', {}).get('datavalue', {})
            value = datavalue.get('value')

            if value is None:
                continue

            result = {'value': value['id'] if isinstance(value, dict) else value}

            # Extract qualifiers
            qualifiers = claim.get('qualifiers', {})
            for qual_prop, qual_claims in qualifiers.items():
                qual_values = []
                for qc in qual_claims:
                    try:
                        qv = qc.get('datavalue', {}).get('value')
                        if qv:
                            qual_values.append(qv['id'] if isinstance(qv, dict) else qv)
                    except (KeyError, TypeError):
                        continue
                if qual_values:
                    result[qual_prop] = qual_values

            # Extract references
            references = []
            for ref in claim.get('references', []):
                ref_data = {}
                for ref_prop, ref_snaks in ref.get('snaks', {}).items():
                    ref_values = []
                    for snak in ref_snaks:
                        try:
                            rv = snak.get('datavalue', {}).get('value')
                            if rv:
                                ref_values.append(rv['id'] if isinstance(rv, dict) else rv)
                        except (KeyError, TypeError):
                            continue
                    if ref_values:
                        ref_data[ref_prop] = ref_values
                if ref_data:
                    references.append(ref_data)

            if references:
                result['references'] = references

            results.append(result)
        except (KeyError, TypeError):
            continue

    return results


# =============================================================================
# ENTITY CLASSIFICATION
# =============================================================================

def classify_entity(instance_qids: Set[str], subclass_qids: Set[str],
                    config: ExtractionConfig) -> str:
    """Classify entity into category based on Q-IDs."""
    all_qids = instance_qids | subclass_qids

    # Check in priority order
    if all_qids & config.drug_qids:
        return 'drug'
    elif all_qids & config.gene_qids:
        return 'gene'
    elif all_qids & config.protein_qids:
        return 'protein'
    elif all_qids & config.disease_qids:
        return 'disease'
    elif all_qids & config.pathway_qids:
        return 'pathway'
    elif all_qids & config.tcm_qids:
        return 'tcm'
    elif all_qids & config.ayurveda_qids:
        return 'ayurveda'
    elif all_qids & config.herb_qids:
        return 'herb'
    elif all_qids & config.supplement_qids:
        return 'supplement'
    elif all_qids & config.chemical_qids:
        return 'chemical'
    else:
        return 'other'


# =============================================================================
# ENTITY EXTRACTION
# =============================================================================

def extract_entity(raw_entity: Dict, config: ExtractionConfig) -> Optional[ExtractedEntity]:
    """Extract and normalize entity from raw Wikidata item."""

    # Skip non-items
    if raw_entity.get('type') != 'item':
        return None

    qid = raw_entity.get('id', '')
    claims = raw_entity.get('claims', {})

    # Get instance-of and subclass-of Q-IDs
    instance_qids = set(get_claim_values(claims, 'P31'))
    subclass_qids = set(get_claim_values(claims, 'P279'))
    all_qids = instance_qids | subclass_qids

    # Check if matches any target Q-ID
    if not (all_qids & config.all_target_qids):
        return None

    # Classify entity
    category = classify_entity(instance_qids, subclass_qids, config)

    # Extract basic info
    labels = raw_entity.get('labels', {})
    descriptions = raw_entity.get('descriptions', {})
    aliases = raw_entity.get('aliases', {})
    sitelinks = raw_entity.get('sitelinks', {})

    label = labels.get('en', {}).get('value', '')
    description = descriptions.get('en', {}).get('value', '')
    alias_list = [a['value'] for a in aliases.get('en', [])]
    wikipedia_url = sitelinks.get('enwiki', {}).get('url')

    # Extract identifiers (single-value properties)
    identifier_props = {
        'P233', 'P234', 'P235', 'P231', 'P274',  # Chemical
        'P662', 'P683', 'P592', 'P715', 'P652', 'P267', 'P486',  # Database IDs
        'P351', 'P353', 'P354', 'P594', 'P639', 'P637',  # Gene
        'P352', 'P591', 'P638',  # Protein
        'P3937', 'P2888', 'P2410',  # Pathway
    }

    identifiers = {}
    for prop in identifier_props:
        values = get_claim_values(claims, prop)
        if values:
            identifiers[prop] = values

    # Extract relationships (with qualifiers/references)
    relationship_props = {
        'P129',   # physically interacts with
        'P2175',  # medical condition treated
        'P769',   # drug interaction
        'P688',   # encodes
        'P702',   # encoded by
        'P680',   # molecular function
        'P681',   # cell component
        'P682',   # biological process
        'P361',   # part of
        'P527',   # has part
        'P2868',  # subject has role
    }

    relationships = {}
    for prop in relationship_props:
        rels = get_claim_with_qualifiers(claims, prop)
        if rels:
            relationships[prop] = rels

    # Filter raw claims to only include relevant properties
    filtered_claims = {
        prop: claims[prop]
        for prop in config.properties
        if prop in claims
    }

    return ExtractedEntity(
        qid=qid,
        label=label,
        description=description,
        aliases=alias_list,
        category=category,
        instance_of=list(instance_qids),
        subclass_of=list(subclass_qids),
        identifiers=identifiers,
        relationships=relationships,
        wikipedia_url=wikipedia_url,
        raw_claims=filtered_claims,
    )


# =============================================================================
# OUTPUT WRITERS
# =============================================================================

class JSONLinesWriter:
    """Write entities to JSON Lines file."""

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.file = open(filepath, 'w', encoding='utf-8')

    def write(self, entity: ExtractedEntity):
        data = asdict(entity)
        self.file.write(json.dumps(data, ensure_ascii=False) + '\n')

    def close(self):
        self.file.close()


class SQLiteWriter:
    """Write entities to SQLite database."""

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.conn = sqlite3.connect(filepath)
        self._create_schema()
        self.batch = []
        self.batch_size = 10000

    def _create_schema(self):
        self.conn.executescript('''
            CREATE TABLE IF NOT EXISTS entities (
                qid TEXT PRIMARY KEY,
                label TEXT,
                description TEXT,
                category TEXT,
                aliases TEXT,
                instance_of TEXT,
                subclass_of TEXT,
                wikipedia_url TEXT,
                data JSON
            );

            CREATE TABLE IF NOT EXISTS identifiers (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                qid TEXT,
                property TEXT,
                value TEXT,
                FOREIGN KEY (qid) REFERENCES entities(qid)
            );

            CREATE TABLE IF NOT EXISTS relationships (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                source_qid TEXT,
                property TEXT,
                target_qid TEXT,
                qualifiers JSON,
                FOREIGN KEY (source_qid) REFERENCES entities(qid)
            );

            CREATE INDEX IF NOT EXISTS idx_entities_category ON entities(category);
            CREATE INDEX IF NOT EXISTS idx_identifiers_property ON identifiers(property);
            CREATE INDEX IF NOT EXISTS idx_identifiers_value ON identifiers(value);
            CREATE INDEX IF NOT EXISTS idx_relationships_source ON relationships(source_qid);
            CREATE INDEX IF NOT EXISTS idx_relationships_target ON relationships(target_qid);
            CREATE INDEX IF NOT EXISTS idx_relationships_property ON relationships(property);
        ''')

    def write(self, entity: ExtractedEntity):
        self.batch.append(entity)
        if len(self.batch) >= self.batch_size:
            self._flush()

    def _flush(self):
        if not self.batch:
            return

        entities_data = []
        identifiers_data = []
        relationships_data = []

        for entity in self.batch:
            entities_data.append((
                entity.qid,
                entity.label,
                entity.description,
                entity.category,
                json.dumps(entity.aliases),
                json.dumps(entity.instance_of),
                json.dumps(entity.subclass_of),
                entity.wikipedia_url,
                json.dumps(entity.raw_claims),
            ))

            for prop, values in entity.identifiers.items():
                for value in values:
                    identifiers_data.append((entity.qid, prop, str(value)))

            for prop, rels in entity.relationships.items():
                for rel in rels:
                    target = rel.get('value', '')
                    qualifiers = {k: v for k, v in rel.items() if k != 'value'}
                    relationships_data.append((
                        entity.qid, prop, target, json.dumps(qualifiers)
                    ))

        self.conn.executemany(
            'INSERT OR REPLACE INTO entities VALUES (?,?,?,?,?,?,?,?,?)',
            entities_data
        )
        self.conn.executemany(
            'INSERT INTO identifiers (qid, property, value) VALUES (?,?,?)',
            identifiers_data
        )
        self.conn.executemany(
            'INSERT INTO relationships (source_qid, property, target_qid, qualifiers) VALUES (?,?,?,?)',
            relationships_data
        )
        self.conn.commit()
        self.batch = []

    def close(self):
        self._flush()
        self.conn.close()


class TSVWriter:
    """Write entities to TSV files (one per category)."""

    def __init__(self, output_dir: str):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.files = {}
        self.writers = {}

    def _get_writer(self, category: str):
        if category not in self.writers:
            filepath = self.output_dir / f'{category}.tsv'
            self.files[category] = open(filepath, 'w', newline='', encoding='utf-8')
            self.writers[category] = csv.writer(self.files[category], delimiter='\t')

            # Write header
            self.writers[category].writerow([
                'qid', 'label', 'description', 'drugbank', 'pubchem', 'chebi',
                'chembl', 'uniprot', 'entrez', 'hgnc_symbol', 'ensembl',
                'smiles', 'atc_code', 'wikipedia_url'
            ])

        return self.writers[category]

    def write(self, entity: ExtractedEntity):
        writer = self._get_writer(entity.category)

        def first_or_empty(values):
            return str(values[0]) if values else ''

        writer.writerow([
            entity.qid,
            entity.label,
            entity.description,
            first_or_empty(entity.identifiers.get('P715', [])),  # DrugBank
            first_or_empty(entity.identifiers.get('P662', [])),  # PubChem
            first_or_empty(entity.identifiers.get('P683', [])),  # ChEBI
            first_or_empty(entity.identifiers.get('P592', [])),  # ChEMBL
            first_or_empty(entity.identifiers.get('P352', [])),  # UniProt
            first_or_empty(entity.identifiers.get('P351', [])),  # Entrez
            first_or_empty(entity.identifiers.get('P353', [])),  # HGNC symbol
            first_or_empty(entity.identifiers.get('P594', [])),  # Ensembl
            first_or_empty(entity.identifiers.get('P233', [])),  # SMILES
            first_or_empty(entity.identifiers.get('P267', [])),  # ATC
            entity.wikipedia_url or '',
        ])

    def close(self):
        for f in self.files.values():
            f.close()


# =============================================================================
# MAIN EXTRACTOR
# =============================================================================

class WikidataExtractor:
    """Main Wikidata extraction orchestrator."""

    def __init__(self, config: ExtractionConfig = None):
        self.config = config or ExtractionConfig()
        self.stats = ExtractionStats()
        self.writers = []

    def add_writer(self, writer):
        """Add output writer."""
        self.writers.append(writer)

    def process(self, input_path: str, progress_interval: int = 100000):
        """Process Wikidata dump."""
        logger.info(f"Processing: {input_path}")
        logger.info(f"Target Q-IDs: {len(self.config.all_target_qids)}")
        logger.info(f"Properties to extract: {len(self.config.properties)}")

        for raw_entity in stream_entities(input_path):
            self.stats.total_processed += 1

            # Progress reporting
            if self.stats.total_processed % progress_interval == 0:
                rate = self.stats.processing_rate()
                logger.info(
                    f"Processed: {self.stats.total_processed:,} | "
                    f"Matched: {self.stats.total_matched:,} | "
                    f"Rate: {rate:.0f}/sec"
                )

            # Extract entity
            entity = extract_entity(raw_entity, self.config)
            if entity is None:
                continue

            self.stats.total_matched += 1
            self.stats.by_category[entity.category] += 1

            # Write to all outputs
            for writer in self.writers:
                writer.write(entity)

        # Close all writers
        for writer in self.writers:
            writer.close()

        # Final statistics
        logger.info("=" * 60)
        logger.info("EXTRACTION COMPLETE")
        logger.info("=" * 60)
        logger.info(f"Total processed: {self.stats.total_processed:,}")
        logger.info(f"Total matched: {self.stats.total_matched:,}")
        logger.info(f"Time elapsed: {self.stats.elapsed_time():.1f} seconds")
        logger.info(f"Processing rate: {self.stats.processing_rate():.0f} entities/sec")
        logger.info("")
        logger.info("By category:")
        for category, count in sorted(self.stats.by_category.items()):
            logger.info(f"  {category}: {count:,}")

        return self.stats


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Extract biomedical entities from Wikidata dump'
    )
    parser.add_argument(
        '--input', '-i', required=True,
        help='Input Wikidata JSON dump path (bz2/gz/json)'
    )
    parser.add_argument(
        '--output-dir', '-o', default='./wikidata_extracted',
        help='Output directory'
    )
    parser.add_argument(
        '--formats', '-f', nargs='+',
        choices=['jsonl', 'sqlite', 'tsv'],
        default=['jsonl', 'sqlite'],
        help='Output formats'
    )
    parser.add_argument(
        '--progress', '-p', type=int, default=100000,
        help='Progress reporting interval'
    )

    args = parser.parse_args()

    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize extractor
    extractor = WikidataExtractor()

    # Add requested writers
    if 'jsonl' in args.formats:
        extractor.add_writer(JSONLinesWriter(output_dir / 'wikidata_biomedical.jsonl'))

    if 'sqlite' in args.formats:
        extractor.add_writer(SQLiteWriter(output_dir / 'wikidata_biomedical.db'))

    if 'tsv' in args.formats:
        extractor.add_writer(TSVWriter(output_dir / 'tsv'))

    # Run extraction
    extractor.process(args.input, progress_interval=args.progress)


if __name__ == '__main__':
    main()
```

### 7.2 ID Mapping Generator

```python
#!/usr/bin/env python3
"""
Generate cross-reference ID mapping files from extracted Wikidata data.

Creates mapping files for:
- Wikidata QID <-> DrugBank
- Wikidata QID <-> PubChem
- Wikidata QID <-> ChEMBL
- Wikidata QID <-> UniProt
- Wikidata QID <-> Entrez Gene
- etc.
"""

import json
import argparse
from pathlib import Path
from collections import defaultdict
from typing import Dict, Set


def create_id_maps(input_jsonl: str, output_dir: str):
    """Create ID mapping files from extracted JSONL."""

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Mapping property codes to names
    id_mappings = {
        'P715': 'drugbank',
        'P662': 'pubchem_cid',
        'P683': 'chebi',
        'P592': 'chembl',
        'P352': 'uniprot',
        'P351': 'entrez_gene',
        'P353': 'hgnc_symbol',
        'P354': 'hgnc_id',
        'P594': 'ensembl_gene',
        'P705': 'ensembl_protein',
        'P231': 'cas',
        'P267': 'atc_code',
        'P652': 'fda_unii',
        'P486': 'mesh',
        'P3937': 'reactome',
        'P2888': 'kegg_pathway',
        'P2410': 'wikipathways',
        'P233': 'smiles',
        'P235': 'inchikey',
    }

    # Initialize mapping dictionaries
    qid_to_id = {name: {} for name in id_mappings.values()}
    id_to_qid = {name: defaultdict(set) for name in id_mappings.values()}

    # Process JSONL
    with open(input_jsonl, 'r') as f:
        for line in f:
            entity = json.loads(line)
            qid = entity['qid']

            for prop, name in id_mappings.items():
                values = entity.get('identifiers', {}).get(prop, [])
                if values:
                    qid_to_id[name][qid] = values
                    for value in values:
                        id_to_qid[name][str(value)].add(qid)

    # Write mapping files
    for name in id_mappings.values():
        # QID -> External ID
        with open(output_path / f'qid_to_{name}.json', 'w') as f:
            json.dump(qid_to_id[name], f, indent=2)

        # External ID -> QID (convert sets to lists for JSON)
        id_to_qid_json = {k: list(v) for k, v in id_to_qid[name].items()}
        with open(output_path / f'{name}_to_qid.json', 'w') as f:
            json.dump(id_to_qid_json, f, indent=2)

    # Create combined mapping file
    combined = {
        'qid_to_ids': {},
        'stats': {}
    }

    with open(input_jsonl, 'r') as f:
        for line in f:
            entity = json.loads(line)
            qid = entity['qid']
            ids = {}

            for prop, name in id_mappings.items():
                values = entity.get('identifiers', {}).get(prop, [])
                if values:
                    ids[name] = values[0] if len(values) == 1 else values

            if ids:
                combined['qid_to_ids'][qid] = ids

    # Stats
    for name in id_mappings.values():
        combined['stats'][name] = len(qid_to_id[name])

    with open(output_path / 'combined_id_map.json', 'w') as f:
        json.dump(combined, f, indent=2)

    print(f"Created ID mapping files in {output_path}")
    print("\nMapping statistics:")
    for name, count in combined['stats'].items():
        print(f"  {name}: {count:,} QIDs")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create ID mapping files')
    parser.add_argument('input', help='Input JSONL file')
    parser.add_argument('--output-dir', '-o', default='./id_maps',
                        help='Output directory')

    args = parser.parse_args()
    create_id_maps(args.input, args.output_dir)
```

### 7.3 Merge with Specialized Databases

```python
#!/usr/bin/env python3
"""
Merge Wikidata extracted data with specialized databases.

Enriches Wikidata entities with:
- DrugBank compound details
- ChEMBL bioactivity data
- UniProt protein details
- BATMAN-TCM traditional medicine data
"""

import json
import sqlite3
import argparse
from pathlib import Path
from typing import Dict, Optional


class DatabaseMerger:
    """Merge Wikidata with specialized databases."""

    def __init__(self, wikidata_db: str):
        self.wd_conn = sqlite3.connect(wikidata_db)
        self.wd_conn.row_factory = sqlite3.Row

    def get_wikidata_by_id(self, property: str, value: str) -> Optional[Dict]:
        """Get Wikidata entity by external identifier."""
        cursor = self.wd_conn.execute('''
            SELECT e.* FROM entities e
            JOIN identifiers i ON e.qid = i.qid
            WHERE i.property = ? AND i.value = ?
        ''', (property, value))

        row = cursor.fetchone()
        if row:
            return dict(row)
        return None

    def merge_drugbank(self, drugbank_file: str, output_file: str):
        """Merge with DrugBank data."""
        merged = []

        # DrugBank format: XML or TSV
        # This is a simplified example - actual implementation
        # would parse DrugBank XML

        with open(drugbank_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue

                drugbank_id = parts[0]

                # Find in Wikidata
                wd_entity = self.get_wikidata_by_id('P715', drugbank_id)

                merged.append({
                    'drugbank_id': drugbank_id,
                    'drugbank_data': parts[1:],
                    'wikidata': wd_entity,
                })

        with open(output_file, 'w') as f:
            json.dump(merged, f, indent=2)

        return len(merged)

    def merge_uniprot(self, uniprot_file: str, output_file: str):
        """Merge with UniProt data."""
        merged = []

        with open(uniprot_file, 'r') as f:
            for line in f:
                if line.startswith('AC'):
                    uniprot_id = line[5:].strip().rstrip(';')
                    wd_entity = self.get_wikidata_by_id('P352', uniprot_id)

                    if wd_entity:
                        merged.append({
                            'uniprot_id': uniprot_id,
                            'wikidata_qid': wd_entity['qid'],
                            'wikidata_label': wd_entity['label'],
                        })

        with open(output_file, 'w') as f:
            json.dump(merged, f, indent=2)

        return len(merged)

    def create_unified_view(self, output_file: str):
        """Create unified view with all identifiers."""
        cursor = self.wd_conn.execute('''
            SELECT
                e.qid,
                e.label,
                e.category,
                GROUP_CONCAT(CASE WHEN i.property = 'P715' THEN i.value END) as drugbank,
                GROUP_CONCAT(CASE WHEN i.property = 'P662' THEN i.value END) as pubchem,
                GROUP_CONCAT(CASE WHEN i.property = 'P683' THEN i.value END) as chebi,
                GROUP_CONCAT(CASE WHEN i.property = 'P592' THEN i.value END) as chembl,
                GROUP_CONCAT(CASE WHEN i.property = 'P352' THEN i.value END) as uniprot,
                GROUP_CONCAT(CASE WHEN i.property = 'P351' THEN i.value END) as entrez,
                GROUP_CONCAT(CASE WHEN i.property = 'P353' THEN i.value END) as hgnc_symbol,
                GROUP_CONCAT(CASE WHEN i.property = 'P594' THEN i.value END) as ensembl
            FROM entities e
            LEFT JOIN identifiers i ON e.qid = i.qid
            GROUP BY e.qid
        ''')

        with open(output_file, 'w') as f:
            f.write('\t'.join([
                'qid', 'label', 'category', 'drugbank', 'pubchem', 'chebi',
                'chembl', 'uniprot', 'entrez', 'hgnc_symbol', 'ensembl'
            ]) + '\n')

            for row in cursor:
                f.write('\t'.join([str(v or '') for v in row]) + '\n')

    def close(self):
        self.wd_conn.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge databases')
    parser.add_argument('wikidata_db', help='Wikidata SQLite database')
    parser.add_argument('--unified', help='Output unified TSV')

    args = parser.parse_args()

    merger = DatabaseMerger(args.wikidata_db)

    if args.unified:
        merger.create_unified_view(args.unified)
        print(f"Created unified view: {args.unified}")

    merger.close()
```

---

## Download

### Official Wikidata Dumps

| Resource | URL | Format | Frequency | Size |
|----------|-----|--------|-----------|------|
| Complete dump | https://dumps.wikimedia.org/wikidatawiki/entities/ | JSON Lines (.jl.bz2) | Weekly | ~90 GB compressed |
| Truthy dump | https://dumps.wikimedia.org/wikidatawiki/entities/ | N-Triples (.nt.bz2) | Weekly | ~20 GB compressed |
| RDF dump | https://dumps.wikimedia.org/wikidatawiki/entities/ | RDF/XML | Weekly | ~40 GB compressed |
| Latest version | https://www.wikidata.org/wiki/Q1 | All formats | Real-time | Varies |

### SPARQL Query Endpoint

```bash
# Query Wikidata SPARQL endpoint (no rate limit for simple queries)
curl -G "https://query.wikidata.org/sparql" \
  --data-urlencode "query=SELECT ?item WHERE { ?item wdt:P31 wd:Q12140 } LIMIT 10" \
  -H "Accept: application/sparql-results+json"
```

### Recommended Download Strategy

```bash
# Download latest truthy dump (smaller, suitable for most use cases)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-truthy.nt.bz2

# For filtered biomedical data, use SPARQL queries instead (more selective)
# See Master SPARQL Queries section for examples
```

---

## Data Format

### JSON Lines Format (.jl.bz2)

Each line is a complete Wikidata item as JSON:

```json
{
  "type": "item",
  "id": "Q12140",
  "labels": {
    "en": { "language": "en", "value": "medication" },
    "de": { "language": "de", "value": "Medikament" }
  },
  "claims": {
    "P31": [
      {
        "rank": "normal",
        "mainsnak": {
          "snaktype": "value",
          "property": "P31",
          "datavalue": { "value": { "entity-type": "item", "numeric-id": 8386 }, "type": "wikibase-entityid" }
        }
      }
    ]
  }
}
```

### Truthy RDF Format (.nt.bz2)

N-Triples RDF format (one statement per line):

```ntriples
<http://www.wikidata.org/entity/Q12140> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://wikiba.se/ontology#Item> .
<http://www.wikidata.org/entity/Q12140> <http://www.w3.org/2000/01/rdf-schema#label> "medication"@en .
<http://www.wikidata.org/entity/Q12140> <http://www.wikidata.org/prop/direct/P31> <http://www.wikidata.org/entity/Q8386> .
```

### SPARQL JSON Results

```json
{
  "head": {
    "vars": ["item", "itemLabel", "drugbank"]
  },
  "results": {
    "bindings": [
      {
        "item": { "type": "uri", "value": "http://www.wikidata.org/entity/Q12140" },
        "itemLabel": { "type": "literal", "value": "medication", "xml:lang": "en" },
        "drugbank": { "type": "literal", "value": "DB00001" }
      }
    ]
  }
}
```

---

## Schema

### Core Wikidata Item Structure

| Field | Type | Description |
|-------|------|-------------|
| `id` | String | Q-ID identifier (Q + digits) |
| `type` | String | Entity type: "item", "property" |
| `labels` | Object | Labels in multiple languages (language→{language, value}) |
| `descriptions` | Object | Descriptions in multiple languages |
| `aliases` | Object | Alternative names in multiple languages |
| `claims` | Object | Claims keyed by P-code; each value is array of statements |
| `sitelinks` | Object | Links to Wikipedia articles (wiki→{site, title}) |

### Claim/Statement Structure

| Field | Type | Description |
|-------|------|-------------|
| `rank` | String | "preferred", "normal", "deprecated" |
| `mainsnak` | Object | Main statement value |
| `qualifiers` | Object | Additional properties (optional) |
| `references` | Array | Citation information (optional) |

### Mainsnak Structure

| Field | Type | Description |
|-------|------|-------------|
| `snaktype` | String | "value", "somevalue", "novalue" |
| `property` | String | P-code property identifier |
| `datavalue` | Object | Typed value (wikibase-entityid, string, quantity, time, etc.) |

### Property (P-code) Relevant to Biomedical Data

| P-code | Label | Datatype | Common Values |
|--------|-------|----------|----------------|
| `P31` | instance of | Item | Q12140 (drug), Q7187 (gene) |
| `P279` | subclass of | Item | Hierarchical parent |
| `P715` | DrugBank ID | String | DB00945 |
| `P661` | ChEMBL ID | String | CHEMBL308 |
| `P2175` | medical condition treated | Item | Q disease_id |
| `P527` | has part | Item | Q compound_id |
| `P352` | UniProt ID | String | P12345 |
| `P486` | MeSH ID | String | D000001 |

---

## Sample Data

### Query: Fetch Medication with External IDs

**SPARQL Query:**
```sparql
SELECT ?drug ?drugLabel ?drugbank ?chebi ?pubchem WHERE {
  ?drug wdt:P31 wd:Q12140 .
  OPTIONAL { ?drug wdt:P715 ?drugbank }
  OPTIONAL { ?drug wdt:P592 ?chebi }
  OPTIONAL { ?drug wdt:P662 ?pubchem }
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en" }
}
LIMIT 5
```

**Sample Results (JSON):**
```json
{
  "head": { "vars": ["drug", "drugLabel", "drugbank", "chebi", "pubchem"] },
  "results": {
    "bindings": [
      {
        "drug": { "type": "uri", "value": "http://www.wikidata.org/entity/Q18216" },
        "drugLabel": { "type": "literal", "value": "aspirin", "xml:lang": "en" },
        "drugbank": { "type": "literal", "value": "DB00945" },
        "chebi": { "type": "literal", "value": "15365" },
        "pubchem": { "type": "literal", "value": "2244" }
      },
      {
        "drug": { "type": "uri", "value": "http://www.wikidata.org/entity/Q12261" },
        "drugLabel": { "type": "literal", "value": "ibuprofen", "xml:lang": "en" },
        "drugbank": { "type": "literal", "value": "DB01050" },
        "chebi": { "type": "literal", "value": "5855" },
        "pubchem": { "type": "literal", "value": "3672" }
      }
    ]
  }
}
```

### Sample Record from Dump (JSON Lines Format)

```json
{
  "type": "item",
  "id": "Q18216",
  "labels": {
    "en": { "language": "en", "value": "aspirin" },
    "de": { "language": "de", "value": "Aspirin" }
  },
  "claims": {
    "P31": [
      {
        "rank": "normal",
        "mainsnak": {
          "snaktype": "value",
          "property": "P31",
          "datavalue": { "value": { "entity-type": "item", "numeric-id": 12140 }, "type": "wikibase-entityid" }
        }
      }
    ],
    "P715": [
      {
        "rank": "normal",
        "mainsnak": {
          "snaktype": "value",
          "property": "P715",
          "datavalue": { "value": "DB00945", "type": "string" }
        }
      }
    ]
  }
}
```

---

## License

### Wikidata License

- **License:** CC0 1.0 Universal (Public Domain Dedication)
- **URL:** https://creativecommons.org/publicdomain/zero/1.0/
- **Summary:** All Wikidata content is in the public domain
- **Attribution:** Optional but appreciated; acknowledge Wikidata as source
- **Restrictions:** None - use for any purpose without limitations
- **Citation:**
  ```
  Wikidata contributors. "Wikidata: The Free Knowledge Base That Anyone Can Edit."
  Retrieved [DATE] from https://www.wikidata.org/
  ```

### Wikipedia Sitelinks License

- **License:** CC-BY-SA 3.0 / CC-BY-SA 4.0 (varies by language)
- **Requirement:** Attribution and ShareAlike when redistributing
- **Note:** Applies only to Wikipedia article text, not Wikidata structured data

### API Terms of Service

- **Rate Limits:** SPARQL endpoint: 60 requests per minute (no official limit for dumps)
- **User-Agent Required:** Identify your application in requests
- **Contact:** wikidata-developers@lists.wikimedia.org for high-volume queries

---

## Data Set Size

### Wikidata Database Statistics

| Metric | Value | Notes |
|--------|-------|-------|
| Total items | ~100 million | As of Jan 2026 |
| Medications (Q12140) | ~45,000 | Primary drug class |
| Pharmaceutical compounds | ~500,000 | All chemical entities |
| Medicinal plants | ~10,000 | Actively used in medicine |
| Properties (P-codes) | ~10,000+ | Schema definitions |

### Storage Estimates

| Dump Type | Compressed | Uncompressed | Notes |
|-----------|-----------|---------------|-------|
| Complete JSON Lines | ~90 GB | ~500 GB | All items and statements |
| Truthy RDF | ~20 GB | ~100 GB | Best-ranked claims only |
| Filtered biomedical | ~500 MB | ~2 GB | Selected Q/P-codes |
| Index database (SQLite) | ~50 GB | ~150 GB | Optimized for queries |

### Processing Time Estimates

| Operation | Time | Hardware |
|-----------|------|----------|
| Download complete dump | 2-4 hours | 100 Mbps connection |
| Parse and index dump | 8-12 hours | 16 CPU cores, 64 GB RAM |
| SPARQL biomedical query | <1 second | Official SPARQL endpoint |
| Custom filtering | 1-2 hours | Extract ~500K biomedical items |

### Last Updated

- **Current Dump:** Weekly (every Monday ~18:00 UTC)
- **SPARQL Endpoint:** Real-time (lag <1 minute)
- **This Reference:** January 2026

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| `Q-ID` | Unique Wikidata item identifier (Q followed by number) | Q12140 (medication) |
| `P-code` | Wikidata property identifier (P followed by number) | P715 (DrugBank ID) |
| `claim` | Statement about an item consisting of property and value | P31=Q12140 |
| `instance_of` | P31 property defining what type an item is | Q12140 instance of Q8386 |
| `subclass_of` | P279 property for hierarchical classification | Q35456 subclass of Q12140 |
| `sitelink` | Link from Wikidata item to Wikipedia article | enwiki:Aspirin |
| `qualifier` | Additional information attached to a claim | start_time, end_time |
| `reference` | Source citation for a claim | stated_in, retrieved |
| `truthy` | Simplified dump containing only best-ranked claims | latest-truthy.nt.bz2 |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| SPARQL | Query language for Wikidata and RDF databases | Wikidata Query Service |
| JSON-LD | JSON format for linked data used in Wikidata dumps | Serialization |
| bz2 | bzip2 compression format for Wikidata dumps | File format |
| jsonl | JSON Lines format with one entity per line | Output format |
| dump | Complete Wikidata database export file | Data extraction |
| streaming_parser | Parser that processes data without loading entirely into memory | Memory efficiency |
| batch_insert | Database operation inserting multiple rows at once | Performance |
| entity | Any item or property in Wikidata | Q-item, P-property |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| QID | Q-Item Identifier | Wikidata entity ID |
| SPARQL | SPARQL Protocol and RDF Query Language | Query interface |
| RDF | Resource Description Framework | Linked data format |
| JSON | JavaScript Object Notation | Data format |
| TSV | Tab-Separated Values | Export format |
| CSV | Comma-Separated Values | Export format |
| SQLite | Self-contained SQL database engine | Local storage |
| API | Application Programming Interface | Programmatic access |
| REST | Representational State Transfer | API architecture |
| CID | Compound Identifier | PubChem ID |
| SMILES | Simplified Molecular Input Line Entry System | Chemical notation |
| InChI | International Chemical Identifier | Chemical notation |
| InChIKey | Hashed InChI for searching | Chemical notation |
| CAS | Chemical Abstracts Service | Registry number |
| UNII | Unique Ingredient Identifier | FDA substance ID |
| ATC | Anatomical Therapeutic Chemical | Drug classification |
| MeSH | Medical Subject Headings | NLM vocabulary |
| HGNC | HUGO Gene Nomenclature Committee | Gene symbols |
| TCM | Traditional Chinese Medicine | Medical system |
| KEGG | Kyoto Encyclopedia of Genes and Genomes | Pathway database |

---

## Appendix A: Quick Reference

### Download Commands

```bash
# Full dump (bzip2 - smaller)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2

# Full dump (gzip - faster decompress)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.gz

# Truthy only (smaller, simpler)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-truthy.nt.bz2

# Fast multi-connection download
aria2c -x 16 -s 16 https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2
```

### Key Q-IDs Quick Reference

| Category | Primary Q-ID | Name |
|----------|--------------|------|
| Drug | Q12140 | medication |
| Chemical | Q11173 | chemical compound |
| Gene | Q7187 | gene |
| Protein | Q8054 | protein |
| Enzyme | Q8047 | enzyme |
| Disease | Q12136 | disease |
| Pathway | Q4915012 | biological pathway |
| TCM | Q891104 | traditional Chinese medicine |
| Herb | Q188617 | medicinal plant |
| Supplement | Q324546 | dietary supplement |

### Key P-codes Quick Reference

| Type | P-code | Name |
|------|--------|------|
| Classification | P31 | instance of |
| Classification | P279 | subclass of |
| Chemical | P233 | SMILES |
| Chemical | P715 | DrugBank ID |
| Chemical | P662 | PubChem CID |
| Gene | P351 | Entrez Gene ID |
| Gene | P353 | HGNC symbol |
| Protein | P352 | UniProt ID |
| Medical | P2175 | treats condition |
| Medical | P129 | interacts with |

---

## Appendix B: Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| Out of memory | Use streaming parser, increase batch size |
| JSON parse error | Skip malformed lines, check encoding |
| Slow processing | Use bz2 (smaller), parallelize with chunks |
| Missing data | Check P31 hierarchy with P279 |
| Duplicate entries | Use INSERT OR REPLACE, deduplicate by QID |

### Performance Tips

1. **Use bz2 format** - Smaller download, acceptable decompression speed
2. **Stream processing** - Never load entire dump into memory
3. **Batch database inserts** - 10,000+ rows per commit
4. **Index after loading** - Create indexes after bulk insert
5. **Filter early** - Check P31/P279 before extracting properties

---

*Document generated: January 2026*
*Last updated: January 18, 2026*
