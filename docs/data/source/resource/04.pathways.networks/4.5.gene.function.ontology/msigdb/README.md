---
id: msigdb
title: "MSigDB - Molecular Signatures Database"
type: source
parent: ../README.md
tier: 1
status: active
category: pathways.networks
subcategory: gene.function.ontology
tags:
  - gene-sets
  - gsea
  - enrichment
  - signatures
  - broad-institute
---

# MSigDB - Molecular Signatures Database

**Category:** [Pathways & Networks](../../_index.md) > [Gene Function & Ontology](../_index.md)

## Overview

MSigDB (Molecular Signatures Database) is a collection of annotated gene sets for use with Gene Set Enrichment Analysis (GSEA) and other enrichment methods. Developed at the Broad Institute, MSigDB provides curated gene sets from pathway databases, published literature, expert knowledge, and expression signatures.

The database contains over 33,000 gene sets organized into 9 major collections covering hallmark signatures, positional gene sets, curated pathways, motif-based sets, computational signatures, GO terms, oncogenic signatures, immunologic signatures, and cell type signatures. MSigDB is the most comprehensive resource for pathway and signature analysis in genomics.

MSigDB integrates gene sets from sources including GO, KEGG, Reactome, WikiPathways, BioCarta, and thousands of published studies.

## Key Statistics

| Metric | Value |
|--------|-------|
| Total Gene Sets | 33,000+ |
| Hallmark Gene Sets | 50 |
| Curated Gene Sets (C2) | 6,300+ |
| GO Gene Sets (C5) | 15,000+ |
| Immunologic Signatures (C7) | 5,200+ |
| Cell Type Signatures (C8) | 700+ |
| Unique Genes | ~40,000 |

## Gene Set Collections

| Collection | ID | Description | Sets |
|------------|----|--------------|----|
| Hallmark | H | Well-defined biological states | 50 |
| Positional | C1 | Chromosome cytogenetic bands | 299 |
| Curated | C2 | Curated pathway databases | 6,300+ |
| Motif | C3 | Transcription factor targets | 3,700+ |
| Computational | C4 | Cancer gene neighborhoods | 850+ |
| Gene Ontology | C5 | GO terms (BP, MF, CC) | 15,000+ |
| Oncogenic | C6 | Oncogenic pathway signatures | 189 |
| Immunologic | C7 | Immunologic signatures | 5,200+ |
| Cell Type | C8 | Cell type markers | 700+ |

## Primary Use Cases

1. **Gene Set Enrichment Analysis (GSEA)** - Identify enriched pathways/signatures
2. **Over-representation analysis (ORA)** - Test gene list overlap
3. **Single-sample enrichment** - Score samples for pathway activity
4. **Biomarker discovery** - Find diagnostic gene signatures
5. **Drug mechanism analysis** - Compare to published drug signatures

## Key Identifiers

| Identifier | Pattern | Example |
|------------|---------|---------|
| Gene Set Name | Collection_Name | HALLMARK_APOPTOSIS |
| Systematic Name | M{#####} | M5930 |
| Gene Symbol | HGNC symbol | TP53 |
| Entrez Gene ID | Numeric | 7157 |

## Hallmark Gene Sets (H Collection)

Representative well-defined biological states:

| Gene Set | Description | Genes |
|----------|-------------|-------|
| HALLMARK_APOPTOSIS | Apoptosis regulation | 161 |
| HALLMARK_HYPOXIA | Hypoxia response | 200 |
| HALLMARK_P53_PATHWAY | p53 signaling | 200 |
| HALLMARK_INFLAMMATORY_RESPONSE | Inflammation | 200 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | EMT | 200 |
| HALLMARK_MYC_TARGETS_V1 | MYC targets | 200 |
| HALLMARK_GLYCOLYSIS | Glycolysis | 200 |
| HALLMARK_DNA_REPAIR | DNA repair | 150 |

## Access Methods

| Method | URL | Notes |
|--------|-----|-------|
| Web Interface | https://www.gsea-msigdb.org/gsea/msigdb/ | Browse/download |
| GSEA Software | https://www.gsea-msigdb.org/gsea/downloads.jsp | Desktop analysis |
| API | https://www.gsea-msigdb.org/gsea/msigdb/api/ | Programmatic access |
| Downloads | https://www.gsea-msigdb.org/gsea/downloads.jsp | GMT files |

### API Examples

```bash
# Search gene sets by keyword
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets?keyword=apoptosis"

# Get gene set details
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/HALLMARK_APOPTOSIS"

# Get genes in a gene set
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets/HALLMARK_APOPTOSIS/genes"

# Search by gene
curl "https://www.gsea-msigdb.org/gsea/msigdb/api/v1/gene_sets?gene=TP53"
```

## Data Formats

| Format | Description | Use Case |
|--------|-------------|----------|
| GMT | Gene Matrix Transposed | GSEA input |
| GMX | Gene Matrix | Alternative format |
| GRP | Gene set list | Simple gene list |
| XML | Detailed metadata | Full annotation |
| JSON | API responses | Programmatic access |

### GMT Format

Tab-delimited: `gene_set_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...`

```
HALLMARK_APOPTOSIS	http://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_APOPTOSIS	CASP3	CASP8	BAX	BCL2	...
```

## Curated Pathway Sources (C2)

| Subcollection | Source | Sets |
|---------------|--------|------|
| CP:BIOCARTA | BioCarta | 292 |
| CP:KEGG | KEGG | 186 |
| CP:PID | NCI-PID | 196 |
| CP:REACTOME | Reactome | 1,615 |
| CP:WIKIPATHWAYS | WikiPathways | 664 |
| CGP | Chemical/Genetic Perturbations | 3,400+ |

## GO Collections (C5)

| Subcollection | Aspect | Sets |
|---------------|--------|------|
| GO:BP | Biological Process | 7,600+ |
| GO:MF | Molecular Function | 1,700+ |
| GO:CC | Cellular Component | 1,000+ |
| HPO | Human Phenotype Ontology | 5,000+ |

## License

| Aspect | Value |
|--------|-------|
| License | CC BY 4.0 |
| Commercial Use | Yes |
| Attribution | Required |
| GSEA Software | Free for non-commercial |

## Cross-References

| Database | Relationship |
|----------|--------------|
| GO | Gene Ontology terms |
| KEGG | Pathway gene sets |
| Reactome | Pathway gene sets |
| WikiPathways | Pathway gene sets |
| PubMed | Literature references |
| GEO | Expression signatures |
| HGNC | Gene symbols |
| Entrez Gene | Gene identifiers |

## GSEA Analysis Steps

1. **Load expression data** - Ranked gene list or expression matrix
2. **Select gene set collection** - Choose relevant MSigDB collection
3. **Run enrichment** - GSEA, ORA, or ssGSEA
4. **Filter results** - FDR q-value < 0.25 typical cutoff
5. **Interpret biology** - Review leading edge genes

## Limitations

- Human-centric; ortholog mapping needed for other species
- Some gene sets derived from older studies may be outdated
- Overlap between gene sets can affect enrichment statistics
- GSEA software has non-commercial license restrictions

## Download Files

| Collection | File |
|------------|------|
| Hallmark | h.all.v2024.1.Hs.symbols.gmt |
| Curated | c2.all.v2024.1.Hs.symbols.gmt |
| GO:BP | c5.go.bp.v2024.1.Hs.symbols.gmt |
| Immunologic | c7.all.v2024.1.Hs.symbols.gmt |
| All | msigdb.v2024.1.Hs.symbols.gmt |

## See Also

- [Gene Ontology](../gene.ontology/_index.md) - Source for GO gene sets
- [Reactome](../../4.1.metabolic.pathways/reactome/_index.md) - Pathway source
- [WikiPathways](../../4.1.metabolic.pathways/wikipathways/_index.md) - Pathway source
