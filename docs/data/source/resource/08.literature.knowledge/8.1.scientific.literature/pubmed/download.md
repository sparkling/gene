---
id: download-pubmed
title: "PubMed Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-23
---

# PubMed Download Instructions

## Quick Start

```bash
# Download PubMed baseline XML (first file)
wget https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed24n0001.xml.gz
```

## Prerequisites

- **wget** or **curl** for downloads
- **Entrez Direct** for E-utilities
- Approximately 100GB-500GB disk space for full baseline

## No Registration Required for Downloads

API key recommended for heavy E-utilities usage (free registration).

## Download Methods

### Method 1: PubMed Baseline (Full Archive)

```bash
# Download complete baseline (annual release)
mkdir pubmed_baseline && cd pubmed_baseline

# Download all baseline files
for i in $(seq -w 1 1219); do
  wget "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/pubmed24n${i}.xml.gz"
done

# Or use rsync (faster)
rsync -av ftp.ncbi.nlm.nih.gov::pubmed/baseline/*.xml.gz ./baseline/
```

### Method 2: Daily Update Files

```bash
# Download daily updates
wget https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/

# Get specific update
wget "https://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/pubmed24n1220.xml.gz"
```

### Method 3: E-utilities API

```bash
# Install Entrez Direct
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# Search PubMed
esearch -db pubmed -query "BRCA1 AND breast cancer" | \
  efetch -format xml > brca1_cancer.xml

# Get specific PMIDs
efetch -db pubmed -id 12345678,23456789 -format xml > specific_articles.xml

# Search with date range
esearch -db pubmed -query "COVID-19" -datetype PDAT -mindate 2023/01/01 -maxdate 2023/12/31 | \
  efetch -format xml > covid_2023.xml

# Get abstract text only
esearch -db pubmed -query "machine learning drug discovery" | \
  efetch -format abstract > ml_drug_abstracts.txt
```

### Method 4: Batch Downloads with E-utilities

```bash
# Large batch download (with API key for higher rate)
export NCBI_API_KEY="your_api_key_here"

# Search and save PMIDs
esearch -db pubmed -query "genomics[MeSH]" | \
  efetch -format uid > genomics_pmids.txt

# Batch fetch by PMID list
cat pmids.txt | \
  xargs -I {} sh -c 'efetch -db pubmed -id {} -format xml >> articles.xml; sleep 0.5'

# Using epost for large batches
cat pmids.txt | epost -db pubmed | efetch -format xml > batch_articles.xml
```

### Method 5: Open Access Subset (Full Text)

```bash
# Download PMC Open Access list
wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_file_list.csv

# Download specific article (full text)
wget "https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/00/00/PMC10000001.tar.gz"

# Bulk download OA articles
rsync -av --include='*.xml' --exclude='*' \
  ftp.ncbi.nlm.nih.gov::pub/pmc/oa_bulk/ ./pmc_oa/
```

### Method 6: Author Manuscript Collection

```bash
# NIH author manuscripts
wget https://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/

# Download by year
rsync -av ftp.ncbi.nlm.nih.gov::pub/pmc/manuscript/2023/ ./manuscripts_2023/
```

## File Inventory

### Baseline Files

| File Pattern | Size | Description |
|--------------|------|-------------|
| pubmed24n*.xml.gz | 10-30 MB each | Baseline XML files |
| Total baseline | ~350 GB | ~36 million citations |

### Update Files

| File Pattern | Frequency | Description |
|--------------|-----------|-------------|
| pubmed24n*.xml.gz | Daily | New/updated records |

### PMC Open Access

| Data Type | Size | Description |
|-----------|------|-------------|
| oa_file_list.csv | ~100 MB | OA article list |
| Full text XML | ~1 TB total | OA full text |

## Post-Download Processing

```bash
# Parse XML with Python
python3 << 'EOF'
import gzip
import xml.etree.ElementTree as ET
import pandas as pd

data = []
with gzip.open('pubmed24n0001.xml.gz', 'rt') as f:
    tree = ET.parse(f)
    root = tree.getroot()

    for article in root.findall('.//PubmedArticle'):
        pmid = article.find('.//PMID').text
        title = article.find('.//ArticleTitle')
        abstract = article.find('.//AbstractText')

        data.append({
            'pmid': pmid,
            'title': title.text if title is not None else None,
            'abstract': abstract.text if abstract is not None else None
        })

df = pd.DataFrame(data)
df.to_csv('pubmed_parsed.tsv', sep='\t', index=False)
print(f"Processed {len(df)} articles")
EOF

# Extract MeSH terms
python3 << 'EOF'
import gzip
import xml.etree.ElementTree as ET

with gzip.open('pubmed24n0001.xml.gz', 'rt') as f:
    tree = ET.parse(f)
    root = tree.getroot()

    with open('mesh_terms.tsv', 'w') as out:
        out.write("pmid\tmesh_term\tmesh_id\n")
        for article in root.findall('.//PubmedArticle'):
            pmid = article.find('.//PMID').text
            for mesh in article.findall('.//MeshHeading'):
                descriptor = mesh.find('DescriptorName')
                if descriptor is not None:
                    term = descriptor.text
                    mesh_id = descriptor.get('UI')
                    out.write(f"{pmid}\t{term}\t{mesh_id}\n")
EOF

# Streaming parse for large files
python3 << 'EOF'
import gzip
from lxml import etree

def parse_pubmed_stream(filename):
    with gzip.open(filename, 'rb') as f:
        context = etree.iterparse(f, events=('end',), tag='PubmedArticle')
        for event, elem in context:
            pmid = elem.find('.//PMID')
            yield {
                'pmid': pmid.text if pmid is not None else None,
                # Add more fields as needed
            }
            elem.clear()

for record in parse_pubmed_stream('pubmed24n0001.xml.gz'):
    print(record['pmid'])
EOF
```

## Verification

```bash
# Check XML structure
zcat pubmed24n0001.xml.gz | head -100

# Count articles
zcat pubmed24n0001.xml.gz | grep -c "<PubmedArticle>"

# Verify E-utilities
esearch -db pubmed -query "test" | head
```

## Update Schedule

| Data Type | Frequency |
|-----------|-----------|
| Baseline release | Annual (December) |
| Update files | Daily |
| API data | Real-time |

## Common Issues

- **Large downloads**: Use rsync for reliability and resume
- **XML parsing**: Use streaming for large files to avoid memory issues
- **Rate limits**: Use API key; limit to 10 req/sec with key, 3/sec without
- **Character encoding**: UTF-8; handle special characters carefully
- **Deleted records**: Update files include deletions; track state

## E-utilities Parameters

| Parameter | Description |
|-----------|-------------|
| db | Database (pubmed, pmc, etc.) |
| query | Search terms |
| rettype | Return type (xml, abstract, etc.) |
| retmax | Maximum records (default 20) |
| datetype | Date field (PDAT, EDAT, etc.) |
| mindate/maxdate | Date range |

## API Key Benefits

| Limit | Without Key | With Key |
|-------|-------------|----------|
| Requests/second | 3 | 10 |
| Batch size | 500 | 500 |
| Email required | Yes | No |

Register at: https://www.ncbi.nlm.nih.gov/account/

## XML Element Reference

```xml
<PubmedArticle>
  <MedlineCitation>
    <PMID>12345678</PMID>
    <Article>
      <ArticleTitle>...</ArticleTitle>
      <Abstract>
        <AbstractText>...</AbstractText>
      </Abstract>
      <AuthorList>...</AuthorList>
    </Article>
    <MeshHeadingList>...</MeshHeadingList>
  </MedlineCitation>
  <PubmedData>
    <ArticleIdList>
      <ArticleId IdType="doi">...</ArticleId>
    </ArticleIdList>
  </PubmedData>
</PubmedArticle>
```

## Related Resources

- [PubMed Central](../pubmed.central/) - Full text articles
- [Europe PMC](../europe.pmc/) - European literature
- [Semantic Scholar](../semantic.scholar/) - AI-enriched literature
