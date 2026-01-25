---
id: download-mesh
title: "Medical Subject Headings (MeSH) Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-23
---

# Medical Subject Headings (MeSH) Download Instructions

## Quick Start

```bash
# Download MeSH descriptors (XML format)
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.xml

# Download supplementary concepts
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/supp2026.xml
```

## Prerequisites

- **wget** or **curl** for downloads
- **XML parser** (lxml, ElementTree)
- Approximately 1GB disk space for full dataset

## No Registration Required

MeSH is public domain and freely available.

## Download Methods

### Method 1: NLM FTP Downloads (XML)

```bash
# MeSH Descriptors (main terms)
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.xml

# Supplementary Concept Records (chemicals, diseases)
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/supp2026.xml

# Qualifiers (subheadings)
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/qual2026.xml

# All files compressed
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/desc2026.gz
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/supp2026.gz
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/xmlmesh/qual2026.gz
```

### Method 2: ASCII Format (Legacy)

```bash
# ASCII descriptors
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2026.bin

# ASCII qualifiers
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/q2026.bin

# ASCII supplementary concepts
wget https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/c2026.bin
```

### Method 3: RDF/Linked Data

```bash
# N-Triples format (large file)
wget https://nlmpubs.nlm.nih.gov/projects/mesh/rdf/mesh.nt.gz

# Download vocabulary
wget https://nlmpubs.nlm.nih.gov/projects/mesh/rdf/vocabulary.nt

# Download 2026 version
wget https://nlmpubs.nlm.nih.gov/projects/mesh/rdf/2026/mesh2026.nt.gz
```

### Method 4: SPARQL Endpoint

```bash
# Query MeSH SPARQL endpoint
curl -X POST "https://id.nlm.nih.gov/mesh/sparql" \
  -H "Content-Type: application/sparql-query" \
  -H "Accept: application/json" \
  -d "SELECT ?d ?name WHERE { ?d a meshv:Descriptor . ?d rdfs:label ?name } LIMIT 100" \
  -o mesh_descriptors.json

# Get disease terms
curl -X POST "https://id.nlm.nih.gov/mesh/sparql" \
  -H "Content-Type: application/sparql-query" \
  -H "Accept: application/json" \
  -d "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
      PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
      SELECT ?d ?name WHERE {
        ?d a meshv:Descriptor .
        ?d rdfs:label ?name .
        ?d meshv:treeNumber ?tree .
        FILTER(STRSTARTS(?tree, 'C'))
      } LIMIT 1000" \
  -o mesh_diseases.json

# Get tree hierarchy for diabetes
curl -X POST "https://id.nlm.nih.gov/mesh/sparql" \
  -H "Content-Type: application/sparql-query" \
  -H "Accept: application/json" \
  -d "PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
      PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
      SELECT ?d ?name ?tree WHERE {
        ?d rdfs:label 'Diabetes Mellitus'@en .
        ?d meshv:treeNumber ?tree
      }" \
  -o mesh_diabetes_tree.json
```

### Method 5: MeSH Browser API

```bash
# Search terms
curl "https://meshb.nlm.nih.gov/api/search/record?searchInField=termDescriptor&sort=&size=25&searchType=exactMatch&searchMethod=SubString&q=diabetes" \
  -o mesh_search.json

# Get descriptor by UI
curl "https://meshb.nlm.nih.gov/api/record/ui/D003920" \
  -o mesh_d003920.json

# Get tree structure
curl "https://meshb.nlm.nih.gov/api/tree/children/C" \
  -o mesh_tree_c.json
```

## File Inventory

### XML Files (Primary)

| File | Size | Description |
|------|------|-------------|
| desc2026.xml | ~300 MB | Main descriptors (~30,000) |
| supp2026.xml | ~500 MB | Supplementary concepts (~280,000) |
| qual2026.xml | ~1 MB | Qualifiers (~80) |

### ASCII Files (Legacy)

| File | Size | Description |
|------|------|-------------|
| d2026.bin | ~100 MB | Descriptors |
| c2026.bin | ~150 MB | Supplementary |
| q2026.bin | ~200 KB | Qualifiers |

### RDF Files

| File | Size | Description |
|------|------|-------------|
| mesh.nt.gz | ~500 MB | Full N-Triples |
| vocabulary.nt | ~50 KB | Vocabulary definitions |

## Post-Download Processing

```bash
# Parse MeSH XML with Python
python3 << 'EOF'
import xml.etree.ElementTree as ET
import gzip

# Parse compressed XML
with gzip.open('desc2026.gz', 'rt', encoding='utf-8') as f:
    tree = ET.parse(f)
    root = tree.getroot()

# Extract descriptors
descriptors = []
for record in root.findall('.//DescriptorRecord'):
    ui = record.find('DescriptorUI').text
    name = record.find('DescriptorName/String').text

    # Get tree numbers
    trees = [tn.text for tn in record.findall('.//TreeNumber')]

    descriptors.append({
        'ui': ui,
        'name': name,
        'trees': ';'.join(trees)
    })

print(f"Total descriptors: {len(descriptors)}")

# Save to TSV
import pandas as pd
df = pd.DataFrame(descriptors)
df.to_csv('mesh_descriptors.tsv', sep='\t', index=False)
EOF

# Extract disease terms (Category C)
python3 << 'EOF'
import xml.etree.ElementTree as ET
import pandas as pd

tree = ET.parse('desc2026.xml')
root = tree.getroot()

diseases = []
for record in root.findall('.//DescriptorRecord'):
    trees = [tn.text for tn in record.findall('.//TreeNumber')]
    # Category C = Diseases
    if any(t.startswith('C') for t in trees):
        ui = record.find('DescriptorUI').text
        name = record.find('DescriptorName/String').text
        diseases.append({
            'ui': ui,
            'name': name,
            'trees': ';'.join(trees)
        })

df = pd.DataFrame(diseases)
df.to_csv('mesh_diseases.tsv', sep='\t', index=False)
print(f"Disease descriptors: {len(df)}")
EOF

# Extract entry terms (synonyms)
python3 << 'EOF'
import xml.etree.ElementTree as ET

tree = ET.parse('desc2026.xml')
root = tree.getroot()

with open('mesh_synonyms.tsv', 'w') as f:
    f.write("descriptor_ui\tdescriptor_name\tentry_term\n")
    for record in root.findall('.//DescriptorRecord'):
        ui = record.find('DescriptorUI').text
        name = record.find('DescriptorName/String').text
        for concept in record.findall('.//Concept'):
            for term in concept.findall('.//Term'):
                term_str = term.find('String').text
                if term_str != name:
                    f.write(f"{ui}\t{name}\t{term_str}\n")
EOF

# Parse ASCII format
python3 << 'EOF'
# MeSH ASCII format is record-based
# Each record starts with *NEWRECORD

records = []
current_record = {}

with open('d2026.bin', 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if line == '*NEWRECORD':
            if current_record:
                records.append(current_record)
            current_record = {}
        elif ' = ' in line:
            key, value = line.split(' = ', 1)
            if key in current_record:
                if isinstance(current_record[key], list):
                    current_record[key].append(value)
                else:
                    current_record[key] = [current_record[key], value]
            else:
                current_record[key] = value

print(f"Records parsed: {len(records)}")
EOF
```

## Verification

```bash
# Check XML structure
head -50 desc2026.xml

# Count descriptors
grep "<DescriptorUI>" desc2026.xml | wc -l

# Count disease terms (tree C)
grep "<TreeNumber>C" desc2026.xml | wc -l

# Check specific term
grep -A 20 "D003920" desc2026.xml | head -30

# Verify RDF format
zcat mesh.nt.gz | head -100
```

## Dataset Versions

| Version | Release Date | Size | Status |
|---------|--------------|------|--------|
| MeSH 2026 | 2025-11-14 | ~800 MB | Current |
| MeSH 2025 | 2024-11 | ~750 MB | Archived |
| Weekly SCR updates | Continuous | Varies | Incremental |

### Version Notes

MeSH 2026 release includes:
- 30,000+ descriptors
- 290,000+ supplementary concept records
- 200,000+ terms including synonyms
- Complete tree hierarchy for 16 categories
- RDF/Linked Data format support

## API Access

| Property | Value |
|----------|-------|
| Base URL | `https://id.nlm.nih.gov/mesh/sparql` |
| Rate Limit | Reasonable use |
| Auth Required | No |
| Documentation | https://meshb.nlm.nih.gov/help |

## Update Schedule

| Release | Frequency |
|---------|-----------|
| Annual release | November (for next year) |
| Weekly updates | Supplementary concepts |
| MeSH on Demand | Continuous indexing |

## Common Issues

- **XML encoding**: Use UTF-8 encoding when parsing
- **Tree numbers**: A descriptor can have multiple tree positions
- **Supplementary vs descriptors**: SCRs are primarily chemicals
- **Entry terms**: Include synonyms, previous names, related terms
- **Major topics**: PubMed uses asterisk for major MeSH terms

## MeSH Tree Structure

| Category | Code | Description |
|----------|------|-------------|
| Anatomy | A | Body structures |
| Organisms | B | Living things |
| Diseases | C | Pathological conditions |
| Chemicals/Drugs | D | Substances |
| Analytical Techniques | E | Methods |
| Psychiatry/Psychology | F | Mental processes |
| Biological Sciences | G | Life science concepts |
| Physical Sciences | H | Non-biological sciences |
| Anthropology/Education | I | Social sciences |
| Technology/Industry | J | Applied sciences |
| Humanities | K | Arts and culture |
| Information Science | L | Communication |
| Named Groups | M | Population categories |
| Healthcare | N | Health services |
| Geographic | Z | Locations |

## Integration Examples

```bash
# Map MeSH to UMLS
# MeSH is included in UMLS Metathesaurus
grep "MSH" MRCONSO.RRF | head -100

# Map MeSH to MONDO
python3 << 'EOF'
import pronto

mondo = pronto.Ontology("mondo.obo")

# Find MeSH cross-references
for term in mondo.terms():
    for xref in term.xrefs:
        if xref.id.startswith("MESH:"):
            print(f"{term.id}\t{term.name}\t{xref.id}")
EOF

# PubMed search with MeSH
# Use MeSH UI in PubMed search
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=D003920[MeSH]&retmax=10" \
  -o pubmed_diabetes_mesh.xml
```

## Related Resources

- [ICD](../icd/) - Clinical classification
- [MONDO](../mondo/) - Disease ontology with MeSH mappings
- [PubMed](../../../../06.literature.knowledge/6.1.scientific.publications/pubmed/) - Literature indexed with MeSH
