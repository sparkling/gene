---
id: integration-xrefs
title: "Integration & Cross-Reference Databases"
type: integration
parent: ../_index.md
last_updated: 2026-01-22
status: migrated
tags: [integration, cross-references, apis]
---

# Integration & Cross-Reference Databases

**Document ID:** 43-86-INTEGRATION-XREFS
**Status:** Final
**Owner:** Data Engineering
**Last Updated:** January 2026
**Version:** 1.0
**Parent:** [../_index.md](../_index.md)

---

## TL;DR

Cross-reference databases are the "glue" connecting genomic, literature, and functional domains. **UniProt ID Mapping** provides the master cross-reference file linking 286 databases with 22-column pre-computed mappings for Gene, RefSeq, PDB, GO, Ensembl, and PubMed. **Gene Ontology** supplies 28 annotation datasets (972K human annotations) linking genes to molecular functions, biological processes, and cellular components. **PMC ID Converter** bridges literature identifiers (PMID/PMCID/DOI) across 200 IDs per request. **Europe PMC Annotations** adds text-mined gene/protein/disease links directly in publications. These form the integration backbone for any multi-domain genetics platform.

---

## Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| **Master ID Hub** | UniProt ID Mapping | 286 databases, pre-computed 22-column file |
| **Functional Annotations** | Gene Ontology GAF 2.2 | Standard format, 28 organism datasets |
| **Literature ID Bridge** | PMC ID Converter | PMID/PMCID/DOI/MID conversion |
| **Text-Mined Links** | Europe PMC Annotations | Pre-extracted gene/protein/disease mentions |
| **Scholarly Metadata** | OpenAlex | CC0, 250M works with PMID/DOI/MAG cross-refs |
| **Knowledge Graph** | Wikidata (see 43-84) | CC0, 286 external identifier properties |

---

## Database Catalog

### 1. UniProt ID Mapping

#### 1.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | UniProt Consortium (SIB, EBI, PIR) |
| **URL** | https://www.uniprot.org/id-mapping |
| **FTP** | ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/ |
| **Total Databases** | 286 cross-referenced databases |
| **Update Frequency** | Every 8 weeks (with UniProt release) |
| **License** | CC BY 4.0 |
| **Commercial Use** | YES (with attribution) |

#### 1.2 ID Mapping File Formats

**idmapping.dat.gz (~15 GB compressed)**

| Column | Content | Example |
|--------|---------|---------|
| 1 | UniProtKB Accession | P04637 |
| 2 | External Database Name | GeneID |
| 3 | External Identifier | 7157 |

**idmapping_selected.tab.gz (~3 GB compressed)**

Pre-computed 22-column table with most-requested mappings:

| Column | Database | Example |
|--------|----------|---------|
| 1 | UniProtKB-AC | P04637 |
| 2 | UniProtKB-ID | P53_HUMAN |
| 3 | GeneID (Entrez) | 7157 |
| 4 | RefSeq | NP_000537.3 |
| 5 | GI | 120407068 |
| 6 | PDB | 1AIE;1C26;1DT7... |
| 7 | GO | GO:0000785;GO:0001046... |
| 8 | UniRef100 | UniRef100_P04637 |
| 9 | UniRef90 | UniRef90_P04637 |
| 10 | UniRef50 | UniRef50_P04637 |
| 11 | UniParc | UPI000002ED67 |
| 12 | PIR | A25224 |
| 13 | NCBI-taxon | 9606 |
| 14 | MIM | 191170 |
| 15 | UniGene | Hs.408312 |
| 16 | PubMed | 11433014;11488916... |
| 17 | EMBL | AB082923;AF052180... |
| 18 | EMBL-CDS | BAC16799.1;AAC12971.1... |
| 19 | Ensembl | ENST00000269305 |
| 20 | Ensembl_TRS | ENST00000269305 |
| 21 | Ensembl_PRO | ENSP00000269305 |
| 22 | Additional PubMed | Additional references |

#### 1.3 Database Categories (286 Total)

**Sequence Databases (6)**
- EMBL, GenBank, DDBJ, RefSeq, CCDS, PIR

**Family & Domain Databases (15)**
- InterPro, Pfam, SMART, PROSITE, Gene3D, HAMAP, CDD, Superfamily, DisProt, MobiDB, PRINTS, PIRSF, PANTHER, TIGRFAMs, NCBIfam

**3D Structure Databases (10)**
- PDB, RCSB-PDB, PDBj, SMR, ModBase, AlphaFoldDB, BMRB, EMDB, PCDDB, SASBDB

**Organism-Specific Databases (30+)**
- MGI, FlyBase, WormBase, SGD, PomBase, TAIR, RGD, ZFIN, Xenbase, Araport, CTD, neXtProt, HPA, dictyBase, EcoGene, CGD, HGNC, VectorBase, EuPathDB, GeneCards, MalaCards

**Protein-Protein Interactions (10)**
- BioGRID, DIP, IntAct, MINT, ComplexPortal, CORUM, STRING, ELM, FunCoup, mentha

**PTM Databases (10)**
- PhosphoSitePlus, CarbonylDB, DEPOD, iPTMnet, SwissPalm, GlyConnect, GlyCosmos, GlyGen, MetOSite, PhosphoNET

**Enzyme & Pathway Databases (15)**
- BioCyc, BRENDA, KEGG, Reactome, SABIO-RK, UniPathway, SignaLink, SIGNOR, PlantReactome, PathBank, CAZy, ESTHER, MoonDB, PeroxiBase, REBASE

**Chemistry Databases (8)**
- ChEMBL, DrugBank, BindingDB, SwissLipids, GuidetoPHARMACOLOGY, DrugCentral, ZINC, PubChem

**Proteomic Databases (10)**
- PRIDE, PeptideAtlas, PaxDb, ProteomicsDB, MassIVE, jPOST, CPTAC, TopDownProteomics, Pumba, ProMEX

**Phylogenomic Databases (8)**
- eggNOG, OMA, OrthoDB, InParanoid, GeneTree, TreeFam, PhylomeDB, PAN-GO

**Genetic Variation Databases (6)**
- dbSNP, ClinGen, GenCC, DisGeNET, DMDM, BioMuta

**Genome Annotation (10)**
- Ensembl, KEGG, GeneID, PATRIC, UCSC, WBParaSite, MANE-Select, Ensembl Genomes, GeneDB, VEuPathDB

#### 1.4 API Access

**REST API Base URL:** `https://rest.uniprot.org/idmapping/`

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/run` | POST | Submit ID mapping job |
| `/status/{jobId}` | GET | Check job status |
| `/results/{jobId}` | GET | Retrieve results |
| `/stream/{jobId}` | GET | Stream large results |

**Request Example:**
```bash
# Submit mapping job
curl -X POST "https://rest.uniprot.org/idmapping/run" \
  -H "Content-Type: application/x-www-form-urlencoded" \
  -d "from=UniProtKB_AC-ID&to=GeneID&ids=P04637,P01308"

# Response: {"jobId": "abc123..."}

# Get results
curl "https://rest.uniprot.org/idmapping/results/abc123"
```

**Rate Limits:**
- Anonymous: 1 request/second
- With email header: 3 requests/second

---

### 2. Gene Ontology (GO)

#### 2.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | Gene Ontology Consortium |
| **URL** | https://geneontology.org/ |
| **Annotations URL** | https://current.geneontology.org/annotations/ |
| **Total Annotations** | 28 organism datasets |
| **Human Annotations** | 972,445 (EBI GOA) |
| **Update Frequency** | Monthly (ontology weekly) |
| **License** | CC BY 4.0 |
| **Commercial Use** | YES (with attribution) |

#### 2.2 Ontology Structure

**Three Root Aspects:**

| Aspect | Code | Description | Terms |
|--------|------|-------------|-------|
| Molecular Function | F | Biochemical activities | ~12,000 |
| Biological Process | P | Larger biological programs | ~30,000 |
| Cellular Component | C | Subcellular locations | ~4,500 |

**Relationship Types:**
- `is_a` - Subclass relationship
- `part_of` - Component relationship
- `regulates` - Regulatory relationship
- `has_part` - Inverse of part_of

#### 2.3 GAF 2.2 Format Specification

**Gene Association File (GAF) 2.2 - 17 Columns:**

| Col | Name | Required | Description | Example |
|-----|------|----------|-------------|---------|
| 1 | DB | Yes | Database source | UniProtKB |
| 2 | DB_Object_ID | Yes | Gene/protein ID | P04637 |
| 3 | DB_Object_Symbol | Yes | Gene symbol | TP53 |
| 4 | Qualifier | Yes | Relation type | enables |
| 5 | GO_ID | Yes | GO term | GO:0003677 |
| 6 | DB:Reference | Yes | Citation | PMID:12345678 |
| 7 | Evidence_Code | Yes | Evidence type | IDA |
| 8 | With/From | Cond | Supporting ID | UniProtKB:Q9Y6K9 |
| 9 | Aspect | Yes | F/P/C | F |
| 10 | DB_Object_Name | No | Full name | Cellular tumor antigen p53 |
| 11 | DB_Object_Synonym | No | Aliases | p53|LFS1 |
| 12 | DB_Object_Type | Yes | Entity type | protein |
| 13 | Taxon | Yes | NCBI taxon | taxon:9606 |
| 14 | Date | Yes | Annotation date | 20240115 |
| 15 | Assigned_By | Yes | Source database | UniProt |
| 16 | Annotation_Extension | No | Context relations | part_of(CL:0000066) |
| 17 | Gene_Product_Form_ID | No | Isoform ID | UniProtKB:P04637-2 |

**Evidence Code Categories:**

| Category | Codes | Description |
|----------|-------|-------------|
| Experimental | EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP | Direct experimental evidence |
| Phylogenetic | IBA, IBD, IKR, IRD | Evolutionary inference |
| Computational | ISS, ISO, ISA, ISM, IGC, RCA | Sequence/structural similarity |
| Author | TAS, NAS | Curator statement |
| Curatorial | IC, ND | Inferred by curator |
| Electronic | IEA | Automated annotation |

**Qualifier Relations (GAF 2.2):**

| Aspect | Qualifier | Meaning |
|--------|-----------|---------|
| F | enables | Gene product performs function |
| F | contributes_to | Partial function contribution |
| P | involved_in | Participates in process |
| P | acts_upstream_of | Upstream regulatory effect |
| P | acts_upstream_of_or_within | Upstream or within |
| C | located_in | Physical location |
| C | is_active_in | Active at location |
| C | part_of | Component of structure |
| C | colocalizes_with | Co-localization |

#### 2.4 Annotation Downloads

**Primary Sources:**

| Organism | File | Annotations | Source |
|----------|------|-------------|--------|
| Human | goa_human.gaf.gz | 972,445 | EBI GOA |
| Mouse | mgi.gaf.gz | 744,389 | MGI |
| Rat | rgd.gaf.gz | 590,055 | RGD |
| Yeast | sgd.gaf.gz | 169,668 | SGD |
| Fly | fb.gaf.gz | 186,975 | FlyBase |
| Worm | wb.gaf.gz | 148,000+ | WormBase |
| Zebrafish | zfin.gaf.gz | 170,000+ | ZFIN |
| Arabidopsis | tair.gaf.gz | 220,000+ | TAIR |

**Download URLs:**
```
https://current.geneontology.org/annotations/goa_human.gaf.gz
https://current.geneontology.org/annotations/mgi.gaf.gz
```

**UniProt Complete Proteomes:**
```
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/
# ~20,000 complete proteomes available
```

#### 2.5 QuickGO API

**Base URL:** `https://www.ebi.ac.uk/QuickGO/services/`

| Endpoint | Purpose | Example |
|----------|---------|---------|
| `/ontology/go/terms/{ids}` | Get GO term details | GO:0003677 |
| `/annotation/search` | Search annotations | gene=TP53 |
| `/ontology/go/terms/{id}/descendants` | Get child terms | GO:0003677 |
| `/ontology/go/terms/{id}/ancestors` | Get parent terms | GO:0003677 |

**Search Parameters:**
- `geneProductId` - UniProt accession
- `geneProductType` - protein, miRNA, etc.
- `taxonId` - NCBI taxonomy ID
- `goId` - GO term ID
- `evidenceCode` - Evidence type filter
- `aspect` - F, P, or C

---

### 3. PMC ID Converter

#### 3.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | NCBI/NLM |
| **URL** | https://pmc.ncbi.nlm.nih.gov/tools/id-converter-api/ |
| **Base API** | https://pmc.ncbi.nlm.nih.gov/tools/idconv/api/v1/articles/ |
| **Supported IDs** | PMID, PMCID, DOI, MID |
| **Max IDs/Request** | 200 |
| **Rate Limit** | 3 requests/second |
| **License** | Public domain |

#### 3.2 Supported ID Types

| ID Type | Format | Example | Description |
|---------|--------|---------|-------------|
| PMID | Numeric | 12345678 | PubMed identifier |
| PMCID | PMC + numeric | PMC1234567 | PubMed Central ID |
| DOI | Standard DOI | 10.1000/xyz123 | Digital Object Identifier |
| MID | Agency prefix | NIHMS123456 | Manuscript ID |
| Versioned PMCID | PMCID.version | PMC1234567.1 | Specific article version |

#### 3.3 API Parameters

| Parameter | Values | Description |
|-----------|--------|-------------|
| `ids` | Comma-separated | Up to 200 IDs |
| `idtype` | pmcid, pmid, mid, doi | Override auto-detection |
| `versions` | yes/no | Include version info |
| `showaiid` | yes/no | Show Article Instance IDs |
| `format` | xml, json, csv, html | Output format |
| `tool` | String | Your application name |
| `email` | Email | Contact email |

#### 3.4 Response Formats

**JSON Response:**
```json
{
  "status": "ok",
  "responseDate": "2026-01-21 10:00:00",
  "request": "ids=12345678;format=json",
  "records": [
    {
      "pmid": "12345678",
      "pmcid": "PMC1234567",
      "doi": "10.1000/xyz123",
      "versions": [
        {"pmcid": "PMC1234567.1", "current": "true"}
      ]
    }
  ]
}
```

**CSV Response Columns:**
- PMID, PMCID, DOI, Version, IsLive, ReleaseDate, ErrorMessage

#### 3.5 Bulk ID Mapping File

**Europe PMC Mapping File:**

| Attribute | Value |
|-----------|-------|
| **URL** | ftp://ftp.ebi.ac.uk/pub/databases/pmc/DOI/PMID_PMCID_DOI.csv.gz |
| **Size** | ~328 MB compressed |
| **Update** | Monthly |
| **Format** | CSV (PMID, PMCID, DOI) |

---

### 4. Europe PMC Annotations API

#### 4.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | EMBL-EBI |
| **URL** | https://europepmc.org/AnnotationsApi |
| **Base API** | https://www.ebi.ac.uk/europepmc/annotations_api/ |
| **Content** | Text-mined entities from 33M+ publications |
| **Format** | JSON-LD (Web Annotations Data Model) |
| **License** | Open (varies by article) |

#### 4.2 Annotation Types

| Type | Source | Description | Cross-References |
|------|--------|-------------|------------------|
| Genes/Proteins | Text mining | Named entity recognition | UniProt, NCBI Gene |
| Chemicals | Text mining | Drug/compound mentions | ChEBI, PubChem |
| Diseases | Text mining | Disease mentions | EFO, DO, MONDO |
| Organisms | Text mining | Species mentions | NCBI Taxonomy |
| GO Terms | Text mining | Functional annotations | Gene Ontology |
| Accession Numbers | Text mining | Database IDs | ENA, UniProt, PDB |

#### 4.3 API Endpoints

**Search Annotations:**
```
GET /annotationsByArticleIds?articleIds={pmid}&format=JSON
```

**Annotation Response Schema:**
```json
{
  "@context": "http://europepmc.org/docs/europepmc-annotation-api-vocab.json",
  "id": "annotation_id",
  "type": "Annotation",
  "creator": "europepmc",
  "body": {
    "type": "SpecificResource",
    "source": "http://purl.uniprot.org/uniprot/P04637"
  },
  "target": {
    "type": "SpecificResource",
    "source": "MED:12345678",
    "selector": {
      "type": "TextQuoteSelector",
      "exact": "p53",
      "prefix": "The tumor suppressor ",
      "suffix": " is mutated in"
    }
  }
}
```

#### 4.4 Cross-Reference URIs

| Entity Type | URI Pattern | Example |
|-------------|-------------|---------|
| UniProt Protein | http://purl.uniprot.org/uniprot/{id} | .../P04637 |
| NCBI Gene | http://identifiers.org/ncbigene/{id} | .../7157 |
| Gene Ontology | http://purl.obolibrary.org/obo/GO_{id} | .../GO_0003677 |
| ChEBI | http://purl.obolibrary.org/obo/CHEBI_{id} | .../CHEBI_15377 |
| EFO Disease | http://www.ebi.ac.uk/efo/EFO_{id} | .../EFO_0000311 |

---

### 5. OpenAlex Cross-References

#### 5.1 Work Object External IDs

| Field | Type | Description |
|-------|------|-------------|
| `ids.openalex` | String | OpenAlex ID (W + number) |
| `ids.doi` | String | DOI with https prefix |
| `ids.pmid` | String | PubMed ID |
| `ids.pmcid` | String | PubMed Central ID |
| `ids.mag` | Integer | Microsoft Academic Graph ID |

**Example:**
```json
{
  "ids": {
    "openalex": "W2741809807",
    "doi": "https://doi.org/10.1038/s41586-020-2649-2",
    "pmid": "32939066",
    "pmcid": "PMC7505446",
    "mag": 2741809807
  }
}
```

#### 5.2 Entity Cross-References

| Entity | External IDs Available |
|--------|------------------------|
| Works | DOI, PMID, PMCID, MAG |
| Authors | ORCID, Scopus, MAG |
| Institutions | ROR, GRID, MAG, Wikidata |
| Sources | ISSN, Wikidata, MAG |
| Concepts | Wikidata, MAG, UMLS |

#### 5.3 Bulk Download

| Resource | URL | Size |
|----------|-----|------|
| Full Snapshot | s3://openalex | ~300 GB |
| Works Only | s3://openalex/data/works | ~250 GB |

**Access:**
```bash
aws s3 sync s3://openalex/data/works --no-sign-request /data/openalex/
```

---

### 6. NCBI Entrez Cross-Links (ELink)

#### 6.1 Overview

| Attribute | Value |
|-----------|-------|
| **Provider** | NCBI |
| **API** | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi |
| **Total Databases** | 38+ Entrez databases |
| **Rate Limit** | 3/sec (anon), 10/sec (API key) |

#### 6.2 Key Database Links

| From | To | Link Name | Description |
|------|----|-----------|-------------|
| PubMed | Gene | pubmed_gene | Cited genes |
| PubMed | Protein | pubmed_protein | Cited proteins |
| PubMed | SNP | pubmed_snp | Associated SNPs |
| PubMed | PMC | pubmed_pmc | Full-text link |
| Gene | PubMed | gene_pubmed | Gene citations |
| Gene | Protein | gene_protein | Gene products |
| Gene | SNP | gene_snp | Gene variants |
| Gene | OMIM | gene_omim | Phenotype links |
| Protein | Gene | protein_gene | Encoding gene |
| Protein | Structure | protein_structure | 3D structures |
| SNP | Gene | snp_gene | SNP location |
| SNP | PubMed | snp_pubmed | SNP citations |
| ClinVar | Gene | clinvar_gene | Clinical variants |
| ClinVar | PubMed | clinvar_pubmed | Evidence |

#### 6.3 ELink Query Examples

**Gene to PubMed:**
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?
  dbfrom=gene&db=pubmed&id=7157&linkname=gene_pubmed
```

**PubMed to Related Databases:**
```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?
  dbfrom=pubmed&id=12345678&cmd=llinks
```

---

## Cross-Reference Integration Matrix

### How Databases Connect

```
                    +-----------+
                    |  PubMed   |
                    |  (PMID)   |
                    +-----+-----+
                          |
         +----------------+----------------+
         |                |                |
         v                v                v
   +-----------+    +-----------+    +-----------+
   |    PMC    |    |  Europe   |    | OpenAlex  |
   |  (PMCID)  |<-->|    PMC    |<-->|   (MAG)   |
   +-----------+    +-----------+    +-----------+
         |                |                |
         |          Text Mining            |
         |                |                |
         v                v                v
   +-----------+    +-----------+    +-----------+
   |   Gene    |    |  UniProt  |<-->|   Gene    |
   |  (Entrez) |<-->|  Protein  |    | Ontology  |
   +-----------+    +-----------+    +-----------+
         |                |                |
         v                v                v
   +-----------+    +-----------+    +-----------+
   |   dbSNP   |    |  PharmGKB |    | Reactome  |
   |   (rsID)  |    |           |    | Pathways  |
   +-----------+    +-----------+    +-----------+
         |                |                |
         +--------->  Wikidata  <----------+
                    (Q-numbers)
```

### ID Mapping Coverage Matrix

| Source | PMID | DOI | Gene ID | UniProt | GO | rsID |
|--------|------|-----|---------|---------|----|----- |
| PubMed | - | Y | Via ELink | Via ELink | N | Via ELink |
| PMC | Y | Y | N | N | N | N |
| OpenAlex | Y | Y | N | N | N | N |
| Europe PMC | Y | Y | Text-mined | Text-mined | Text-mined | N |
| UniProt | Y | N | Y | - | Y | N |
| Gene Ontology | Y | N | Y | Y | - | N |
| Wikidata | Y | Y | Y | Y | Partial | Y |

---

## Integration Recommendations

### Phase 1: Core ID Infrastructure (Week 1-2)

| Task | Source | Records | Priority |
|------|--------|---------|----------|
| Download UniProt idmapping_selected.tab | FTP | Full DB | Critical |
| Download Europe PMC ID mappings | FTP | ~30M | Critical |
| Build PMID-DOI-PMCID lookup table | PMC/Europe PMC | ~30M | Critical |
| Index by all ID types | Local | - | Critical |

### Phase 2: Functional Annotations (Week 3-4)

| Task | Source | Records | Priority |
|------|--------|---------|----------|
| Download human GO annotations | GO Consortium | 972K | High |
| Download mouse GO annotations | MGI | 744K | High |
| Map UniProt to GO terms | GAF files | 1.7M+ | High |
| Build GO term hierarchy | OBO file | 46K terms | High |

### Phase 3: Text-Mined Links (Week 5-6)

| Task | Source | Records | Priority |
|------|--------|---------|----------|
| Query Europe PMC annotations | API | On-demand | Medium |
| Cache gene/protein mentions | Local | Variable | Medium |
| Link papers to genes | Computed | Variable | Medium |

### Phase 4: Knowledge Graph Enrichment (Week 7-8)

| Task | Source | Records | Priority |
|------|--------|---------|----------|
| Extract Wikidata cross-refs | SPARQL | ~1M | Medium |
| Build ELink relationship cache | NCBI API | On-demand | Medium |
| Integrate OpenAlex metadata | S3/API | 250M | Low |

---

## Technical Specifications

### RuVector Cross-Reference Schema

```typescript
// ID Mapping Collection
const idMappingCollection = {
  name: 'id_mappings',
  properties: {
    primary_id: { type: 'string', indexed: true },
    primary_type: { type: 'string', indexed: true },  // pmid, uniprot, gene_id, etc.
    mappings: {
      type: 'object',
      properties: {
        pmid: { type: 'string[]', indexed: true },
        pmcid: { type: 'string[]', indexed: true },
        doi: { type: 'string[]', indexed: true },
        uniprot: { type: 'string[]', indexed: true },
        gene_id: { type: 'string[]', indexed: true },
        ensembl: { type: 'string[]', indexed: true },
        hgnc: { type: 'string[]', indexed: true },
        refseq: { type: 'string[]', indexed: true },
        go_terms: { type: 'string[]', indexed: true },
        pdb: { type: 'string[]' }
      }
    }
  }
};

// GO Annotation Collection
const goAnnotationCollection = {
  name: 'go_annotations',
  properties: {
    uniprot_id: { type: 'string', indexed: true },
    gene_symbol: { type: 'string', indexed: true },
    go_id: { type: 'string', indexed: true },
    aspect: { type: 'string', indexed: true },  // F, P, C
    qualifier: { type: 'string' },
    evidence_code: { type: 'string', indexed: true },
    reference: { type: 'string', indexed: true },
    taxon: { type: 'number', indexed: true },
    date: { type: 'date' }
  }
};

// Graph Relationships
// (:Gene)-[:ANNOTATED_WITH]->(:GOTerm)
// (:Protein)-[:ENCODED_BY]->(:Gene)
// (:Article)-[:MENTIONS]->(:Gene)
// (:Article)-[:MENTIONS]->(:Protein)
```

### Storage Estimates

| Dataset | Raw Size | Indexed Size | Notes |
|---------|----------|--------------|-------|
| UniProt ID Mapping (selected) | 3 GB | 5 GB | Indexed by all ID types |
| PMID-DOI-PMCID mappings | 500 MB | 1 GB | ~30M records |
| Human GO annotations | 100 MB | 300 MB | 972K annotations |
| GO Ontology (OBO) | 50 MB | 150 MB | 46K terms + relations |
| **Total** | **~4 GB** | **~7 GB** | |

---

## Dependencies

### Upstream Dependencies

| Dependency | Purpose | Risk if Unavailable |
|------------|---------|---------------------|
| UniProt FTP | ID mapping master file | HIGH - no equivalent |
| GO Consortium | Functional annotations | HIGH - critical for function |
| PMC ID Converter | Literature ID bridge | MEDIUM - Europe PMC backup |
| Europe PMC | Text-mined annotations | LOW - can text-mine locally |
| OpenAlex | Scholarly metadata | LOW - PubMed alternative |

### Downstream Dependents

| Dependent | Usage |
|-----------|-------|
| SNP Evidence System | Gene-to-literature links |
| RAG System | Paper retrieval and context |
| Knowledge Graph | Entity relationships |
| Search | Cross-database ID resolution |
| API | ID conversion endpoints |

---

## Cost Summary

### One-Time Setup

| Item | Cost |
|------|------|
| UniProt download | $0 (FTP) |
| GO download | $0 (FTP) |
| Europe PMC mappings | $0 (FTP) |
| Processing compute | ~$20 |
| **Total** | **~$20** |

### Ongoing Monthly

| Item | Cost |
|------|------|
| Storage (self-hosted) | $0 (included in VPS) |
| API calls | $0 (bulk download) |
| Update processing | $0 (cron) |
| **Total** | **~$0** |

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| UniProt release delay | Low | Medium | Use previous release |
| GO format change | Very Low | Medium | Version-specific parsers |
| API rate limiting | Medium | Low | Bulk downloads preferred |
| ID mapping gaps | Medium | Medium | Multiple source fallbacks |
| Cross-ref staleness | Low | Low | Regular sync schedule |

---

## Glossary

| Term | Definition | Example |
|------|------------|---------|
| Cross-reference (xref) | A link between identifiers in different databases that refer to the same biological entity | UniProt AC P04637 maps to Entrez Gene 7157 |
| ID Mapping | The process of converting identifiers from one database to equivalent identifiers in another | Converting rs IDs to gene symbols |
| Bulk Download | Retrieving large datasets via FTP/HTTP rather than individual API calls | Downloading UniProt idmapping.dat.gz |

### Domain-Specific Terms

| Term | Definition | Related To |
|------|------------|------------|
| UniProt AC | Accession number uniquely identifying a protein in UniProt database | P04637, Q9Y243 |
| Entrez Gene ID | NCBI's unique numeric identifier for genes | 7157 (TP53) |
| HGNC ID | HUGO Gene Nomenclature Committee's approved gene symbol | HGNC:11998 (TP53) |
| Ensembl ID | Stable identifier from Ensembl genome database | ENSG00000141510 |
| GAF 2.2 | Gene Association File format version 2.2 for GO annotations | goa_human.gaf |
| PMC ID | PubMed Central identifier for full-text articles | PMC3531190 |
| GO Term | Gene Ontology identifier describing molecular function, biological process, or cellular component | GO:0006915 |

### Acronyms

| Acronym | Expansion | Notes |
|---------|-----------|-------|
| API | Application Programming Interface | REST, SPARQL endpoints |
| FTP | File Transfer Protocol | Used for bulk downloads |
| GO | Gene Ontology | Controlled vocabulary for gene products |
| HGNC | HUGO Gene Nomenclature Committee | Official human gene naming authority |
| PMC | PubMed Central | Free full-text archive |
| VPS | Virtual Private Server | Self-hosted deployment |
| CC BY | Creative Commons Attribution | License type |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Data Engineering | Initial document |
