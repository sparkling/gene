# Wikidata, Wikipedia, and DBpedia Bulk Download Guide

## Overview

This document provides comprehensive guidance for bulk downloading and processing data from Wikidata, Wikipedia, and DBpedia for pharmaceutical, gene, and protein research applications.

**Last Updated:** January 2026

---

## 1. Wikidata Dumps

### 1.1 Full JSON Dump Location

**Primary URL:** https://dumps.wikimedia.org/wikidatawiki/entities/

### 1.2 Available File Formats and Sizes (as of January 2026)

| File | Format | Compressed Size | Notes |
|------|--------|-----------------|-------|
| `latest-all.json.bz2` | JSON (bzip2) | ~99.6 GB | **Recommended** - Best compression ratio |
| `latest-all.json.gz` | JSON (gzip) | ~151.2 GB | Faster decompression than bz2 |
| `latest-all.nt.bz2` | N-Triples (bzip2) | ~189.8 GB | RDF format |
| `latest-all.nt.gz` | N-Triples (gzip) | ~245.8 GB | RDF format |
| `latest-all.ttl.bz2` | Turtle (bzip2) | ~121.5 GB | RDF format, more readable |
| `latest-all.ttl.gz` | Turtle (gzip) | ~149.1 GB | RDF format |

**Uncompressed Size Estimate:** ~1 TB for JSON format

### 1.3 Truthy vs Full Statements

| Dump Type | File | Size | Use Case |
|-----------|------|------|----------|
| **Full Statements** | `latest-all.json.bz2` | ~99.6 GB | Complete data with all qualifiers, references, and deprecated statements |
| **Truthy Statements** | `latest-truthy.nt.bz2` | ~42.2 GB | Only "best" (truthy) values - simpler, smaller, faster to process |

**Truthy dumps** contain only statements with the highest rank (preferred > normal > deprecated). Use truthy dumps when:
- You need only current/best values
- Storage or processing time is limited
- You don't need statement provenance or qualifiers

### 1.4 Incremental vs Full Dumps

| Type | Frequency | Location | Use Case |
|------|-----------|----------|----------|
| **Full Dumps** | Weekly | `latest-all.json.bz2` | Initial data load, complete refresh |
| **Daily Snapshots** | Daily | `/wikidatawiki/entities/YYYYMMDD/` | Track specific dated versions |
| **Recent Changes API** | Real-time | MediaWiki API | Keep local copy up-to-date |

### 1.5 Lexeme Dumps (For Linguistic Data)

| File | Size | Content |
|------|------|---------|
| `latest-lexemes.json.bz2` | ~421 MB | Words, forms, and senses |
| `latest-lexemes.nt.bz2` | ~1.06 GB | RDF N-Triples format |

### 1.6 Key Entity Classes for Drug/Gene/Protein Filtering

| Entity Type | Q-ID | Description |
|-------------|------|-------------|
| **Medication/Drug** | Q12140 | Pharmaceutical drugs |
| **Chemical Compound** | Q11173 | Chemical substances |
| **Gene** | Q7187 | Basic unit of heredity |
| **Protein** | Q8054 | Biomolecules |
| **Protein-Coding Gene** | Q20747295 | Genes that encode proteins |
| **Disease** | Q12136 | Medical conditions |
| **Biological Process** | Q2996394 | GO biological processes |

### 1.7 Key Properties for Drug-Gene Relationships

| Property | P-ID | Description |
|----------|------|-------------|
| **Instance of** | P31 | Type classification |
| **Subclass of** | P279 | Class hierarchy |
| **Encodes** | P688 | Gene to protein relationship |
| **Encoded by** | P702 | Protein to gene relationship |
| **Ortholog** | P684 | Orthologous genes |
| **Gene Ontology Molecular Function** | P680 | GO:MF annotations |
| **Gene Ontology Cell Component** | P681 | GO:CC annotations |
| **Gene Ontology Biological Process** | P682 | GO:BP annotations |
| **Treats** | P2175 | Drug treats disease |
| **Drug Interaction** | P769 | Drug-drug interactions |
| **Active Ingredient** | P3781 | Drug composition |

---

## 2. Wikipedia Dumps

### 2.1 Database Dumps Location

**Primary URL:** https://dumps.wikimedia.org/

**English Wikipedia:** https://dumps.wikimedia.org/enwiki/latest/

### 2.2 Key Files

| File | Size (enwiki) | Content |
|------|---------------|---------|
| `enwiki-latest-pages-articles.xml.bz2` | ~24.8 GB | Current article revisions only |
| `enwiki-latest-pages-articles-multistream.xml.bz2` | ~25.9 GB | Multistream format (27 parts) |
| `enwiki-latest-all-titles-in-ns0.gz` | ~106.8 MB | Main namespace titles only |
| `enwiki-latest-all-titles.gz` | ~375.5 MB | All namespace titles |

**Uncompressed Size (enwiki):** ~105+ GB for pages-articles

### 2.3 Size Estimates by Language

| Language | Wiki Code | Compressed Size | Articles |
|----------|-----------|-----------------|----------|
| English | enwiki | ~24.8 GB | ~6.8M |
| German | dewiki | ~7.5 GB | ~2.8M |
| French | frwiki | ~6.2 GB | ~2.5M |
| Spanish | eswiki | ~4.8 GB | ~1.9M |
| Japanese | jawiki | ~4.2 GB | ~1.4M |
| Russian | ruwiki | ~5.1 GB | ~1.9M |
| Chinese | zhwiki | ~3.8 GB | ~1.3M |
| Italian | itwiki | ~4.5 GB | ~1.8M |

### 2.4 Dump Schedule

- **Frequency:** Twice monthly (1st and 20th)
- **Availability:** 7-10 days after scheduled date
- **Retention:** Multiple versions available (latest/ symlink to most recent)

### 2.5 Parsing Tools

#### WikiExtractor

**Installation:**
```bash
pip install wikiextractor
```

**Basic Usage:**
```bash
wikiextractor enwiki-latest-pages-articles.xml.bz2 -o output/ --json
```

**Key Options:**
| Option | Description |
|--------|-------------|
| `-o OUTPUT` | Output directory |
| `--json` | JSON output format |
| `--no-templates` | Skip template expansion (faster) |
| `-b n[KMG]` | Max file size (default 1M) |
| `--processes N` | Parallel processes |

**Output Format:**
```xml
<doc id="123" url="https://..." title="Article Title">
    Extracted plain text content...
</doc>
```

#### mwparserfromhell

**Installation:**
```bash
pip install mwparserfromhell
```

**Usage:**
```python
import mwparserfromhell

wikicode = mwparserfromhell.parse("'''Bold''' and [[link]]")
templates = wikicode.filter_templates()
links = wikicode.filter_wikilinks()
text = wikicode.strip_code()
```

**Key Features:**
- Parse MediaWiki wikitext directly
- Extract templates, links, headings
- Compatible with Python 3.9+
- Fast C tokenizer extension

---

## 3. DBpedia Downloads

### 3.1 DBpedia Databus

**Primary URL:** https://databus.dbpedia.org/

DBpedia Databus provides structured access to DBpedia datasets with:
- Monthly release cycles
- Modular dataset organization
- SPARQL-based data selection
- Version control and metadata

### 3.2 Relevant Datasets

#### Mapping-Based Objects
**URL:** `https://databus.dbpedia.org/dbpedia/mappings/mappingbased-objects`

| Aspect | Description |
|--------|-------------|
| **Content** | High-quality RDF statements with IRI object values |
| **Source** | Wikipedia infoboxes via mapping extraction |
| **Quality** | URI canonicalization + type consistency filtering |
| **License** | CC-BY 3.0 |

#### Infobox Properties (Generic)
**URL:** `https://databus.dbpedia.org/dbpedia/generic/infobox-properties`

| Aspect | Description |
|--------|-------------|
| **Content** | All infobox data (less curated) |
| **Coverage** | Best fact coverage |
| **Trade-off** | Less consistency than mapping-based |
| **License** | CC-BY 3.0 |

### 3.3 Available Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| **Turtle** | `.ttl` | Human-readable RDF |
| **N-Triples** | `.nt` | Line-based RDF (streaming-friendly) |
| **N-Quads** | `.nq` | N-Triples with named graphs |

**Note:** Current releases prioritize TTL/TQL formats; NT/NQ may be deprecated.

### 3.4 Download Methods

**Using DBpedia Databus Client:**
```bash
# Install
pip install databus-client

# Query and download
databus-client query "SELECT ?file WHERE { ... }" --output ./data/
```

**Direct HTTP Download:**
```bash
wget https://databus.dbpedia.org/dbpedia/mappings/mappingbased-objects/2022.12.01/mappingbased-objects_lang=en.ttl.bz2
```

---

## 4. Processing Strategies

### 4.1 Streaming Parsers for Large Files

#### JSON Streaming with ijson
```python
import ijson
import bz2

with bz2.open('latest-all.json.bz2', 'rt') as f:
    for item in ijson.items(f, 'item'):
        process_entity(item)
```

#### Memory-Efficient Line-by-Line JSON
```python
import gzip
import json

def stream_wikidata_dump(filepath):
    """Stream Wikidata JSON dump line by line."""
    opener = gzip.open if filepath.endswith('.gz') else bz2.open

    with opener(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            # Skip array brackets
            if line.strip() in ['[', ']']:
                continue
            # Remove trailing comma
            if line.rstrip().endswith(','):
                line = line.rstrip()[:-1]
            try:
                yield json.loads(line)
            except json.JSONDecodeError:
                continue
```

#### XML Streaming with iterparse
```python
from lxml import etree

def stream_wikipedia_dump(filepath):
    """Stream Wikipedia XML dump with memory efficiency."""
    context = etree.iterparse(
        filepath,
        events=('end',),
        tag='{http://www.mediawiki.org/xml/export-0.10/}page'
    )

    for event, elem in context:
        title = elem.findtext('{http://www.mediawiki.org/xml/export-0.10/}title')
        text = elem.findtext('.//{http://www.mediawiki.org/xml/export-0.10/}text')

        yield {'title': title, 'text': text}

        # Critical: Clear element to free memory
        elem.clear()
        # Also clear parent references
        while elem.getprevious() is not None:
            del elem.getparent()[0]
```

### 4.2 Filtering Strategies

#### Filter by Q-ID (Instance Of Property)

```python
def is_pharmaceutical_entity(entity):
    """Check if entity is a drug, gene, or protein."""
    target_types = {
        'Q12140',   # medication
        'Q11173',   # chemical compound
        'Q7187',    # gene
        'Q8054',    # protein
        'Q12136',   # disease
    }

    claims = entity.get('claims', {})
    instance_of = claims.get('P31', [])

    for claim in instance_of:
        try:
            qid = claim['mainsnak']['datavalue']['value']['id']
            if qid in target_types:
                return True
        except (KeyError, TypeError):
            continue

    return False
```

#### Filter by Property Presence

```python
def has_drug_gene_properties(entity):
    """Check if entity has drug-gene relationship properties."""
    relevant_properties = {
        'P688',   # encodes
        'P702',   # encoded by
        'P2175',  # treats
        'P769',   # drug interaction
        'P680',   # GO molecular function
        'P681',   # GO cellular component
        'P682',   # GO biological process
    }

    claims = entity.get('claims', {})
    return bool(relevant_properties & set(claims.keys()))
```

### 4.3 Chunked Processing

```python
from qwikidata.json_dump import WikidataJsonDump

# Create manageable chunks
wjd = WikidataJsonDump('latest-all.json.bz2')

# Split into 100,000 entity chunks
chunk_files = wjd.create_chunks(
    num_lines_per_chunk=100_000,
    output_dir='./chunks/'
)

# Process chunks in parallel
from multiprocessing import Pool

def process_chunk(chunk_file):
    results = []
    for entity in WikidataJsonDump(chunk_file):
        if is_pharmaceutical_entity(entity):
            results.append(extract_data(entity))
    return results

with Pool(8) as pool:
    all_results = pool.map(process_chunk, chunk_files)
```

### 4.4 Resource Requirements

| Dataset | Disk (Compressed) | Disk (Working) | RAM (Streaming) | RAM (Full Load) |
|---------|-------------------|----------------|-----------------|-----------------|
| Wikidata Full JSON | 100 GB | 200 GB | 4-8 GB | N/A (too large) |
| Wikidata Truthy | 42 GB | 100 GB | 2-4 GB | N/A |
| Wikipedia EN | 25 GB | 150 GB | 4-8 GB | N/A |
| DBpedia Mappings | 2-5 GB | 10 GB | 1-2 GB | 16-32 GB |

**Recommended Minimum Setup:**
- **CPU:** 8+ cores for parallel processing
- **RAM:** 16 GB (32 GB preferred)
- **Disk:** 500 GB SSD (1 TB recommended)
- **Network:** 100+ Mbps for downloads

---

## 5. Recommended Pipeline

```
Download  ->  Validate  ->  Stream Parse  ->  Filter  ->  Extract  ->  Transform  ->  Store
   |            |              |               |            |             |           |
   v            v              v               v            v             v           v
 wget/       checksum      ijson/         Q-ID/P-ID     Entity        Normalize    SQLite/
 aria2c     validation    iterparse       matching      fields        schema      PostgreSQL
```

### 5.1 Pipeline Implementation

```python
#!/usr/bin/env python3
"""
Wikidata Bulk Processing Pipeline for Drug-Gene Data
"""

import bz2
import json
import sqlite3
from pathlib import Path
from typing import Iterator, Dict, Any
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor

@dataclass
class DrugGeneRelation:
    drug_qid: str
    drug_label: str
    gene_qid: str
    gene_label: str
    relation_type: str
    references: list

class WikidataProcessor:
    """Process Wikidata dumps for pharmaceutical data."""

    DRUG_TYPES = {'Q12140', 'Q11173', 'Q35456'}  # medication, chemical compound, pharmaceutical drug
    GENE_TYPES = {'Q7187', 'Q20747295'}  # gene, protein-coding gene
    PROTEIN_TYPES = {'Q8054'}  # protein

    def __init__(self, dump_path: str, db_path: str):
        self.dump_path = Path(dump_path)
        self.db_path = Path(db_path)
        self._init_database()

    def _init_database(self):
        """Initialize SQLite database with schema."""
        conn = sqlite3.connect(self.db_path)
        conn.executescript('''
            CREATE TABLE IF NOT EXISTS entities (
                qid TEXT PRIMARY KEY,
                entity_type TEXT,
                label_en TEXT,
                description_en TEXT,
                aliases TEXT,
                data JSON
            );

            CREATE TABLE IF NOT EXISTS drug_gene_relations (
                id INTEGER PRIMARY KEY,
                drug_qid TEXT,
                gene_qid TEXT,
                relation_property TEXT,
                qualifiers JSON,
                references JSON,
                FOREIGN KEY (drug_qid) REFERENCES entities(qid),
                FOREIGN KEY (gene_qid) REFERENCES entities(qid)
            );

            CREATE INDEX IF NOT EXISTS idx_entity_type ON entities(entity_type);
            CREATE INDEX IF NOT EXISTS idx_drug_qid ON drug_gene_relations(drug_qid);
            CREATE INDEX IF NOT EXISTS idx_gene_qid ON drug_gene_relations(gene_qid);
        ''')
        conn.close()

    def stream_entities(self) -> Iterator[Dict[str, Any]]:
        """Stream entities from compressed dump."""
        opener = bz2.open if str(self.dump_path).endswith('.bz2') else open

        with opener(self.dump_path, 'rt', encoding='utf-8') as f:
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
                    print(f"Warning: JSON decode error at line {line_num}: {e}")
                    continue

    def get_entity_type(self, entity: Dict) -> str:
        """Determine entity type from P31 (instance of) claims."""
        claims = entity.get('claims', {})
        p31_claims = claims.get('P31', [])

        for claim in p31_claims:
            try:
                qid = claim['mainsnak']['datavalue']['value']['id']
                if qid in self.DRUG_TYPES:
                    return 'drug'
                elif qid in self.GENE_TYPES:
                    return 'gene'
                elif qid in self.PROTEIN_TYPES:
                    return 'protein'
            except (KeyError, TypeError):
                continue

        return None

    def extract_label(self, entity: Dict, lang: str = 'en') -> str:
        """Extract label in specified language."""
        labels = entity.get('labels', {})
        label_data = labels.get(lang, {})
        return label_data.get('value', '')

    def extract_description(self, entity: Dict, lang: str = 'en') -> str:
        """Extract description in specified language."""
        descriptions = entity.get('descriptions', {})
        desc_data = descriptions.get(lang, {})
        return desc_data.get('value', '')

    def extract_drug_gene_relations(self, entity: Dict) -> list:
        """Extract drug-gene relationship claims."""
        relations = []
        claims = entity.get('claims', {})

        # Key properties for drug-gene relationships
        relation_properties = {
            'P688': 'encodes',           # gene encodes protein
            'P702': 'encoded_by',        # protein encoded by gene
            'P2175': 'treats',           # drug treats disease
            'P769': 'drug_interaction',  # drug-drug interaction
            'P3489': 'mechanism_of_action',
            'P129': 'physically_interacts_with',
        }

        for prop_id, relation_name in relation_properties.items():
            prop_claims = claims.get(prop_id, [])
            for claim in prop_claims:
                try:
                    target_qid = claim['mainsnak']['datavalue']['value']['id']

                    # Extract qualifiers
                    qualifiers = {}
                    for qual_prop, qual_claims in claim.get('qualifiers', {}).items():
                        qualifiers[qual_prop] = [
                            q.get('datavalue', {}).get('value')
                            for q in qual_claims
                        ]

                    # Extract references
                    references = []
                    for ref in claim.get('references', []):
                        ref_data = {}
                        for ref_prop, ref_snaks in ref.get('snaks', {}).items():
                            ref_data[ref_prop] = [
                                s.get('datavalue', {}).get('value')
                                for s in ref_snaks
                            ]
                        references.append(ref_data)

                    relations.append({
                        'source_qid': entity['id'],
                        'target_qid': target_qid,
                        'property': prop_id,
                        'relation_name': relation_name,
                        'qualifiers': qualifiers,
                        'references': references,
                    })
                except (KeyError, TypeError):
                    continue

        return relations

    def process(self, batch_size: int = 10000):
        """Process entire dump and store relevant data."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        entity_batch = []
        relation_batch = []
        processed = 0
        relevant = 0

        print("Starting Wikidata dump processing...")

        for entity in self.stream_entities():
            processed += 1

            if processed % 100000 == 0:
                print(f"Processed: {processed:,} | Relevant: {relevant:,}")

            entity_type = self.get_entity_type(entity)
            if not entity_type:
                continue

            relevant += 1

            # Extract entity data
            entity_data = {
                'qid': entity['id'],
                'entity_type': entity_type,
                'label_en': self.extract_label(entity),
                'description_en': self.extract_description(entity),
                'aliases': json.dumps([
                    a['value'] for a in entity.get('aliases', {}).get('en', [])
                ]),
                'data': json.dumps({
                    'claims': entity.get('claims', {}),
                    'sitelinks': entity.get('sitelinks', {}),
                })
            }
            entity_batch.append(entity_data)

            # Extract relations
            relations = self.extract_drug_gene_relations(entity)
            relation_batch.extend(relations)

            # Batch insert
            if len(entity_batch) >= batch_size:
                self._insert_entities(cursor, entity_batch)
                self._insert_relations(cursor, relation_batch)
                conn.commit()
                entity_batch = []
                relation_batch = []

        # Insert remaining
        if entity_batch:
            self._insert_entities(cursor, entity_batch)
            self._insert_relations(cursor, relation_batch)
            conn.commit()

        conn.close()
        print(f"Complete! Processed: {processed:,} | Stored: {relevant:,}")

    def _insert_entities(self, cursor, entities):
        """Batch insert entities."""
        cursor.executemany('''
            INSERT OR REPLACE INTO entities
            (qid, entity_type, label_en, description_en, aliases, data)
            VALUES (:qid, :entity_type, :label_en, :description_en, :aliases, :data)
        ''', entities)

    def _insert_relations(self, cursor, relations):
        """Batch insert relations."""
        cursor.executemany('''
            INSERT INTO drug_gene_relations
            (drug_qid, gene_qid, relation_property, qualifiers, references)
            VALUES (:source_qid, :target_qid, :property,
                    :qualifiers, :references)
        ''', [
            {**r, 'qualifiers': json.dumps(r['qualifiers']),
             'references': json.dumps(r['references'])}
            for r in relations
        ])


if __name__ == '__main__':
    processor = WikidataProcessor(
        dump_path='latest-all.json.bz2',
        db_path='wikidata_pharma.db'
    )
    processor.process()
```

---

## 6. Code Examples

### 6.1 Download with Progress

```python
#!/usr/bin/env python3
"""Download Wikidata dump with progress tracking."""

import requests
from tqdm import tqdm
from pathlib import Path

def download_wikidata_dump(
    url: str = "https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2",
    output_path: str = "latest-all.json.bz2"
):
    """Download Wikidata dump with progress bar."""

    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))

    output = Path(output_path)

    with open(output, 'wb') as f:
        with tqdm(total=total_size, unit='B', unit_scale=True, desc=output.name) as pbar:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))

    print(f"Downloaded: {output} ({output.stat().st_size / 1e9:.2f} GB)")

# Alternative: Use aria2c for faster multi-connection download
# aria2c -x 16 -s 16 https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2
```

### 6.2 Filter for Pharmaceutical Items (P31 = Q12140)

```python
#!/usr/bin/env python3
"""Filter Wikidata dump for pharmaceutical items."""

import bz2
import json
from typing import Iterator, Dict

def stream_pharmaceutical_entities(dump_path: str) -> Iterator[Dict]:
    """
    Stream only pharmaceutical-related entities from Wikidata dump.

    Filters for:
    - Q12140: medication/pharmaceutical drug
    - Q11173: chemical compound
    - Q35456: pharmaceutical drug
    - Q7187: gene
    - Q8054: protein
    """

    PHARMA_TYPES = {
        'Q12140',   # medication
        'Q11173',   # chemical compound
        'Q35456',   # pharmaceutical drug
        'Q7187',    # gene
        'Q8054',    # protein
        'Q12136',   # disease
        'Q20747295',  # protein-coding gene
    }

    with bz2.open(dump_path, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line in ['[', ']', '']:
                continue
            if line.endswith(','):
                line = line[:-1]

            try:
                entity = json.loads(line)
            except json.JSONDecodeError:
                continue

            # Check P31 (instance of) claims
            claims = entity.get('claims', {})
            p31_claims = claims.get('P31', [])

            for claim in p31_claims:
                try:
                    type_qid = claim['mainsnak']['datavalue']['value']['id']
                    if type_qid in PHARMA_TYPES:
                        yield entity
                        break
                except (KeyError, TypeError):
                    continue


def main():
    """Process and count pharmaceutical entities."""
    counts = {'medication': 0, 'compound': 0, 'gene': 0, 'protein': 0, 'disease': 0}

    for entity in stream_pharmaceutical_entities('latest-all.json.bz2'):
        qid = entity['id']
        label = entity.get('labels', {}).get('en', {}).get('value', 'N/A')

        # Categorize
        claims = entity.get('claims', {})
        p31_claims = claims.get('P31', [])

        for claim in p31_claims:
            try:
                type_qid = claim['mainsnak']['datavalue']['value']['id']
                if type_qid == 'Q12140':
                    counts['medication'] += 1
                elif type_qid == 'Q11173':
                    counts['compound'] += 1
                elif type_qid == 'Q7187':
                    counts['gene'] += 1
                elif type_qid == 'Q8054':
                    counts['protein'] += 1
                elif type_qid == 'Q12136':
                    counts['disease'] += 1
            except (KeyError, TypeError):
                continue

        # Sample output
        if sum(counts.values()) <= 10:
            print(f"{qid}: {label}")

    print(f"\nCounts: {counts}")


if __name__ == '__main__':
    main()
```

### 6.3 Extract Drug-Gene Relationships

```python
#!/usr/bin/env python3
"""Extract drug-gene relationships from Wikidata."""

import bz2
import json
import csv
from typing import Iterator, Dict, List, Tuple
from dataclasses import dataclass, asdict

@dataclass
class DrugGeneRelation:
    drug_qid: str
    drug_label: str
    gene_qid: str
    gene_label: str
    protein_qid: str
    protein_label: str
    relation_type: str
    mechanism: str
    references: str

class DrugGeneExtractor:
    """Extract drug-gene-protein relationships."""

    # Relevant properties
    ENCODES = 'P688'           # gene encodes protein
    ENCODED_BY = 'P702'        # protein encoded by gene
    TREATS = 'P2175'           # treats disease
    MECHANISM = 'P3489'        # mechanism of action (target)
    INTERACTS = 'P129'         # physically interacts with
    INHIBITS = 'P3364'         # inhibits
    ACTIVATES = 'P3771'        # activates

    def __init__(self, dump_path: str):
        self.dump_path = dump_path
        self.entities = {}  # Cache for QID -> label lookup

    def build_entity_index(self, limit: int = None) -> Dict[str, Dict]:
        """Build index of relevant entities (drugs, genes, proteins)."""
        print("Building entity index...")

        RELEVANT_TYPES = {'Q12140', 'Q11173', 'Q7187', 'Q8054', 'Q20747295'}

        count = 0
        with bz2.open(self.dump_path, 'rt', encoding='utf-8') as f:
            for line in f:
                if limit and count >= limit:
                    break

                line = line.strip()
                if line in ['[', ']', '']:
                    continue
                if line.endswith(','):
                    line = line[:-1]

                try:
                    entity = json.loads(line)
                except json.JSONDecodeError:
                    continue

                # Check if relevant type
                claims = entity.get('claims', {})
                is_relevant = False
                entity_type = None

                for claim in claims.get('P31', []):
                    try:
                        qid = claim['mainsnak']['datavalue']['value']['id']
                        if qid in RELEVANT_TYPES:
                            is_relevant = True
                            if qid in {'Q12140', 'Q11173'}:
                                entity_type = 'drug'
                            elif qid in {'Q7187', 'Q20747295'}:
                                entity_type = 'gene'
                            elif qid == 'Q8054':
                                entity_type = 'protein'
                            break
                    except (KeyError, TypeError):
                        continue

                if is_relevant:
                    self.entities[entity['id']] = {
                        'label': entity.get('labels', {}).get('en', {}).get('value', ''),
                        'type': entity_type,
                        'claims': claims,
                    }
                    count += 1

                    if count % 50000 == 0:
                        print(f"Indexed {count:,} entities")

        print(f"Total indexed: {count:,}")
        return self.entities

    def extract_relations(self) -> List[DrugGeneRelation]:
        """Extract drug-gene relationships."""
        relations = []

        for qid, data in self.entities.items():
            if data['type'] != 'drug':
                continue

            drug_label = data['label']
            claims = data['claims']

            # Check mechanism of action (P3489) - links to proteins/genes
            for claim in claims.get(self.MECHANISM, []):
                try:
                    target_qid = claim['mainsnak']['datavalue']['value']['id']
                    if target_qid in self.entities:
                        target = self.entities[target_qid]

                        # Get references
                        refs = self._extract_references(claim)

                        if target['type'] == 'protein':
                            # Find the encoding gene
                            gene_qid, gene_label = self._find_encoding_gene(target_qid)

                            relations.append(DrugGeneRelation(
                                drug_qid=qid,
                                drug_label=drug_label,
                                gene_qid=gene_qid or '',
                                gene_label=gene_label or '',
                                protein_qid=target_qid,
                                protein_label=target['label'],
                                relation_type='mechanism_of_action',
                                mechanism='target',
                                references=refs,
                            ))
                        elif target['type'] == 'gene':
                            relations.append(DrugGeneRelation(
                                drug_qid=qid,
                                drug_label=drug_label,
                                gene_qid=target_qid,
                                gene_label=target['label'],
                                protein_qid='',
                                protein_label='',
                                relation_type='mechanism_of_action',
                                mechanism='gene_target',
                                references=refs,
                            ))
                except (KeyError, TypeError):
                    continue

            # Check physical interactions (P129)
            for claim in claims.get(self.INTERACTS, []):
                try:
                    target_qid = claim['mainsnak']['datavalue']['value']['id']
                    if target_qid in self.entities:
                        target = self.entities[target_qid]
                        refs = self._extract_references(claim)

                        if target['type'] in ('protein', 'gene'):
                            relations.append(DrugGeneRelation(
                                drug_qid=qid,
                                drug_label=drug_label,
                                gene_qid=target_qid if target['type'] == 'gene' else '',
                                gene_label=target['label'] if target['type'] == 'gene' else '',
                                protein_qid=target_qid if target['type'] == 'protein' else '',
                                protein_label=target['label'] if target['type'] == 'protein' else '',
                                relation_type='physical_interaction',
                                mechanism='binding',
                                references=refs,
                            ))
                except (KeyError, TypeError):
                    continue

        return relations

    def _find_encoding_gene(self, protein_qid: str) -> Tuple[str, str]:
        """Find the gene that encodes a protein."""
        if protein_qid not in self.entities:
            return None, None

        protein_data = self.entities[protein_qid]
        for claim in protein_data['claims'].get(self.ENCODED_BY, []):
            try:
                gene_qid = claim['mainsnak']['datavalue']['value']['id']
                if gene_qid in self.entities:
                    return gene_qid, self.entities[gene_qid]['label']
            except (KeyError, TypeError):
                continue

        return None, None

    def _extract_references(self, claim: Dict) -> str:
        """Extract reference information from claim."""
        refs = []
        for ref in claim.get('references', []):
            snaks = ref.get('snaks', {})
            # P248 = stated in (source)
            for snak in snaks.get('P248', []):
                try:
                    source_qid = snak['datavalue']['value']['id']
                    refs.append(source_qid)
                except (KeyError, TypeError):
                    pass
        return '|'.join(refs)

    def save_to_csv(self, relations: List[DrugGeneRelation], output_path: str):
        """Save relations to CSV."""
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=DrugGeneRelation.__dataclass_fields__.keys())
            writer.writeheader()
            for rel in relations:
                writer.writerow(asdict(rel))
        print(f"Saved {len(relations)} relations to {output_path}")


def main():
    extractor = DrugGeneExtractor('latest-all.json.bz2')

    # Build index (optionally limit for testing)
    extractor.build_entity_index(limit=None)  # Set limit=100000 for testing

    # Extract relations
    relations = extractor.extract_relations()

    # Save results
    extractor.save_to_csv(relations, 'drug_gene_relations.csv')

    # Summary
    print(f"\nExtracted {len(relations)} drug-gene relationships")

    # Sample output
    print("\nSample relations:")
    for rel in relations[:5]:
        print(f"  {rel.drug_label} -> {rel.protein_label or rel.gene_label} ({rel.relation_type})")


if __name__ == '__main__':
    main()
```

### 6.4 Using qwikidata Library

```python
#!/usr/bin/env python3
"""Process Wikidata dump using qwikidata library."""

from qwikidata.json_dump import WikidataJsonDump
from qwikidata.entity import WikidataItem

def process_with_qwikidata(dump_path: str):
    """
    Process Wikidata dump using qwikidata library.

    Installation: pip install qwikidata
    """
    wjd = WikidataJsonDump(dump_path)

    # Target types
    DRUG_QID = 'Q12140'
    GENE_QID = 'Q7187'
    PROTEIN_QID = 'Q8054'

    drugs = []
    genes = []
    proteins = []

    for entity_dict in wjd:
        if entity_dict['type'] != 'item':
            continue

        entity = WikidataItem(entity_dict)

        # Check instance of (P31)
        claim_group = entity.get_claim_group('P31')

        for claim in claim_group:
            target = claim.mainsnak.datavalue
            if target is None:
                continue

            target_qid = target.value['id']

            if target_qid == DRUG_QID:
                drugs.append({
                    'qid': entity.entity_id,
                    'label': entity.get_label('en'),
                    'description': entity.get_description('en'),
                })
            elif target_qid == GENE_QID:
                genes.append({
                    'qid': entity.entity_id,
                    'label': entity.get_label('en'),
                })
            elif target_qid == PROTEIN_QID:
                proteins.append({
                    'qid': entity.entity_id,
                    'label': entity.get_label('en'),
                })

        # Progress
        total = len(drugs) + len(genes) + len(proteins)
        if total % 10000 == 0 and total > 0:
            print(f"Found: {len(drugs)} drugs, {len(genes)} genes, {len(proteins)} proteins")

    return drugs, genes, proteins


# Create chunks for parallel processing
def create_chunks(dump_path: str, output_dir: str, lines_per_chunk: int = 100000):
    """Split dump into manageable chunks."""
    wjd = WikidataJsonDump(dump_path)
    chunk_files = wjd.create_chunks(
        num_lines_per_chunk=lines_per_chunk,
        output_dir=output_dir
    )
    print(f"Created {len(chunk_files)} chunks")
    return chunk_files
```

---

## 7. SPARQL Queries for Validation

Use these queries on https://query.wikidata.org/ to validate your extracted data.

### 7.1 Count Pharmaceutical Drugs

```sparql
SELECT (COUNT(?drug) AS ?count) WHERE {
  ?drug wdt:P31 wd:Q12140 .
}
```

### 7.2 Drugs with Gene Targets

```sparql
SELECT ?drug ?drugLabel ?gene ?geneLabel WHERE {
  ?drug wdt:P31 wd:Q12140 .
  ?drug wdt:P129 ?target .
  ?target wdt:P702 ?gene .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

### 7.3 Drug-Disease-Gene Network

```sparql
SELECT DISTINCT ?drug ?drugLabel ?disease ?diseaseLabel ?gene ?geneLabel WHERE {
  ?drug wdt:P31 wd:Q12140 .
  ?drug wdt:P2175 ?disease .
  ?disease wdt:P2293 ?gene .
  SERVICE wikibase:label { bd:serviceParam wikibase:language "en". }
}
LIMIT 100
```

---

## 8. References and Resources

### Official Documentation
- [Wikidata Database Download](https://www.wikidata.org/wiki/Wikidata:Database_download)
- [Wikimedia Downloads](https://dumps.wikimedia.org/)
- [DBpedia Databus](https://databus.dbpedia.org/)

### Tools and Libraries
- [qwikidata](https://github.com/kensho-technologies/qwikidata) - Python Wikidata tools
- [WikiExtractor](https://github.com/attardi/wikiextractor) - Wikipedia text extraction
- [mwparserfromhell](https://github.com/earwig/mwparserfromhell) - MediaWiki parser

### Tutorials and Guides
- [Wikidata Streaming Dump Processing](https://kokes.github.io/blog/2017/11/12/wikidata-streaming-dump.html)
- [simple-wikidata-db](https://github.com/neelguha/simple-wikidata-db) - Preprocessing scripts
- [Wikidata SPARQL Tutorial](https://www.wikidata.org/wiki/Wikidata:SPARQL_tutorial)

### Biological Data in Wikidata
- [Gene Wiki Data Model](https://www.wikidata.org/wiki/Wikidata:WikiProject_Molecular_biology)
- [Drug SPARQL Examples](https://github.com/sebotic/SPARQL/blob/master/drugs_sparql_examples.md)
- [Wikidata as Semantic Framework for Gene Wiki](https://pmc.ncbi.nlm.nih.gov/articles/PMC4795929/)

---

## Appendix A: Quick Download Commands

```bash
# Wikidata JSON dump (bz2 - smaller)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2

# Wikidata JSON dump (gz - faster to decompress)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.gz

# Wikidata truthy statements (smaller, simpler)
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-truthy.nt.bz2

# English Wikipedia articles
wget https://dumps.wikimedia.org/enwiki/latest/enwiki-latest-pages-articles.xml.bz2

# Multi-connection download with aria2c (much faster)
aria2c -x 16 -s 16 https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2

# Verify checksum
wget https://dumps.wikimedia.org/wikidatawiki/entities/latest-all.json.bz2.md5
md5sum -c latest-all.json.bz2.md5
```

## Appendix B: Disk Space Planning

| Step | Space Required | Notes |
|------|----------------|-------|
| Download (compressed) | 100-150 GB | Wikidata JSON |
| Working space | 200-300 GB | Temporary files, chunks |
| Extracted database | 50-100 GB | SQLite/PostgreSQL |
| Wikipedia EN | 25 GB compressed | 150 GB working |
| DBpedia mappings | 5-10 GB | All languages |
| **Total recommended** | **500 GB - 1 TB** | SSD recommended |
