# Research Papers ETL Pipeline Design

**Last Updated:** January 2026
**Target:** 1-10M genetics-relevant papers from PubMed/PMC
**Vector Database:** RuVector (384d embeddings, tiered compression)
**Relational Database:** Supabase PostgreSQL
**Budget Constraint:** Bootstrapped startup, cost-sensitive

---

## Executive Summary

| Metric | Value |
|--------|-------|
| **Source** | PubMed/MEDLINE + PMC Open Access |
| **Total Papers in PubMed** | 36M+ citations |
| **Genetics-Relevant Estimate** | 3-5M papers |
| **Filtered Target** | 1-3M high-quality papers |
| **Storage (RuVector)** | ~400MB-1.2GB (with tiered compression) |
| **Processing Time** | 2-4 days (initial), 1-2 hours (daily updates) |
| **Annual Cost** | ~$0 (all open access, self-hosted processing) |

---

## 1. Download Strategy

### 1.1 PubMed Data Architecture

```
PubMed Data Distribution
├── Baseline Files (Annual)
│   ├── ftp://ftp.ncbi.nlm.nih.gov/pubmed/baseline/
│   ├── 1,219 compressed XML files (~100 GB total)
│   ├── Each file: ~30,000 citations
│   └── Released annually (December)
│
├── Daily Update Files
│   ├── ftp://ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/
│   ├── New citations, corrections, retractions
│   └── ~1,000-2,000 new papers/day
│
└── PMC Open Access Subset
    ├── ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/
    ├── Full-text articles under CC licenses
    └── ~4M articles with full text
```

### 1.2 Download Strategy

#### Primary: Baseline + Incremental Updates

```python
# scripts/etl/extract/pubmed_downloader.py

import os
import ftplib
import gzip
import hashlib
from pathlib import Path
from datetime import datetime, timedelta
import time
from typing import Generator
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class PubMedDownloader:
    """Download PubMed baseline and update files."""

    FTP_HOST = "ftp.ncbi.nlm.nih.gov"
    BASELINE_PATH = "/pubmed/baseline/"
    UPDATE_PATH = "/pubmed/updatefiles/"

    # Rate limiting: NCBI allows 3 requests/second without API key
    REQUESTS_PER_SECOND = 2  # Conservative

    def __init__(self, output_dir: Path, api_key: str = None):
        self.output_dir = output_dir
        self.api_key = api_key  # Increases limit to 10 req/sec
        self.rate_limit = 10 if api_key else 2

        # Create directory structure
        (output_dir / "baseline").mkdir(parents=True, exist_ok=True)
        (output_dir / "updates").mkdir(parents=True, exist_ok=True)
        (output_dir / "checksums").mkdir(parents=True, exist_ok=True)

    def download_baseline(self, max_files: int = None) -> Generator[Path, None, None]:
        """
        Download PubMed baseline files.

        For cost-sensitive operation, we filter during download:
        - Only download files that haven't been downloaded
        - Verify checksums to avoid re-downloading corrupt files
        """
        logger.info("Connecting to NCBI FTP server...")

        with ftplib.FTP(self.FTP_HOST) as ftp:
            ftp.login()  # Anonymous login
            ftp.cwd(self.BASELINE_PATH)

            # Get list of files
            files = []
            ftp.retrlines('NLST', files.append)

            # Filter to only .xml.gz files
            xml_files = sorted([f for f in files if f.endswith('.xml.gz')])

            if max_files:
                xml_files = xml_files[:max_files]

            logger.info(f"Found {len(xml_files)} baseline files to process")

            for i, filename in enumerate(xml_files):
                output_path = self.output_dir / "baseline" / filename

                # Skip if already downloaded and verified
                if output_path.exists() and self._verify_checksum(ftp, filename, output_path):
                    logger.info(f"[{i+1}/{len(xml_files)}] Skipping {filename} (already downloaded)")
                    yield output_path
                    continue

                # Download file
                logger.info(f"[{i+1}/{len(xml_files)}] Downloading {filename}...")

                with open(output_path, 'wb') as f:
                    ftp.retrbinary(f'RETR {filename}', f.write)

                # Rate limiting
                time.sleep(1 / self.rate_limit)

                yield output_path

    def download_updates(self, since_date: datetime = None) -> Generator[Path, None, None]:
        """
        Download daily update files since a given date.

        Default: last 7 days of updates
        """
        if since_date is None:
            since_date = datetime.now() - timedelta(days=7)

        logger.info(f"Downloading updates since {since_date.date()}")

        with ftplib.FTP(self.FTP_HOST) as ftp:
            ftp.login()
            ftp.cwd(self.UPDATE_PATH)

            files = []
            ftp.retrlines('NLST', files.append)

            # Filter to update files (pubmedYYn####.xml.gz format)
            update_files = sorted([f for f in files if f.endswith('.xml.gz')])

            for filename in update_files:
                output_path = self.output_dir / "updates" / filename

                # Check file modification time
                mod_time = ftp.sendcmd(f'MDTM {filename}')[4:]
                file_date = datetime.strptime(mod_time, '%Y%m%d%H%M%S')

                if file_date < since_date:
                    continue

                if output_path.exists():
                    logger.info(f"Skipping {filename} (already downloaded)")
                    yield output_path
                    continue

                logger.info(f"Downloading update {filename}...")

                with open(output_path, 'wb') as f:
                    ftp.retrbinary(f'RETR {filename}', f.write)

                time.sleep(1 / self.rate_limit)

                yield output_path

    def _verify_checksum(self, ftp: ftplib.FTP, filename: str, local_path: Path) -> bool:
        """Verify file checksum against NCBI's MD5 file."""
        md5_filename = filename + ".md5"

        try:
            # Download MD5 checksum
            md5_content = []
            ftp.retrlines(f'RETR {md5_filename}', md5_content.append)
            expected_md5 = md5_content[0].split()[1]

            # Calculate local file MD5
            with open(local_path, 'rb') as f:
                actual_md5 = hashlib.md5(f.read()).hexdigest()

            return expected_md5 == actual_md5

        except Exception as e:
            logger.warning(f"Could not verify checksum for {filename}: {e}")
            return False


class PubMedUpdateTracker:
    """Track downloaded files and last update time."""

    def __init__(self, state_file: Path):
        self.state_file = state_file
        self.state = self._load_state()

    def _load_state(self) -> dict:
        if self.state_file.exists():
            import json
            with open(self.state_file) as f:
                return json.load(f)
        return {
            "baseline_version": None,
            "last_baseline_file": None,
            "last_update_date": None,
            "downloaded_files": []
        }

    def save_state(self):
        import json
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2, default=str)

    def mark_downloaded(self, filename: str):
        if filename not in self.state["downloaded_files"]:
            self.state["downloaded_files"].append(filename)
        self.save_state()

    def is_downloaded(self, filename: str) -> bool:
        return filename in self.state["downloaded_files"]
```

### 1.3 Rate Limiting and API Etiquette

| Parameter | Without API Key | With API Key |
|-----------|-----------------|--------------|
| **E-utilities rate limit** | 3 requests/second | 10 requests/second |
| **FTP download rate** | No official limit | No official limit |
| **Recommended FTP rate** | 2-3 files/second | 5-10 files/second |
| **Daily E-utilities quota** | 100,000 requests | Unlimited |

**Best Practices:**

```python
# Rate limiting configuration
RATE_CONFIG = {
    "ftp_delay_seconds": 0.5,      # 2 files/second
    "api_delay_seconds": 0.4,       # 2.5 requests/second (conservative)
    "batch_size": 100,              # E-utilities batch size
    "retry_attempts": 3,
    "retry_backoff": 2.0,           # Exponential backoff multiplier
    "timeout_seconds": 30,
}

# Required headers for API requests
HEADERS = {
    "tool": "gene-kb-etl",
    "email": "admin@yourdomain.com",  # Required by NCBI
}
```

### 1.4 Storage for Raw Downloads

```
/home/ubuntu/src/gene/data/
├── raw/
│   └── ncbi/
│       └── pubmed/
│           ├── baseline/           # ~100 GB for full baseline
│           │   ├── pubmed24n0001.xml.gz
│           │   ├── pubmed24n0002.xml.gz
│           │   └── ...
│           ├── updates/            # ~5 GB (rolling 30 days)
│           │   ├── pubmed24n1201.xml.gz
│           │   └── ...
│           ├── state.json          # Download tracking
│           └── checksums/
│               └── verified.json
```

**Storage Optimization:**

```python
# Delete baseline files after processing to save space
# Keep only the extracted/filtered data

def cleanup_after_processing(file_path: Path, keep_days: int = 7):
    """Remove processed files older than keep_days."""
    file_age = datetime.now() - datetime.fromtimestamp(file_path.stat().st_mtime)
    if file_age.days > keep_days:
        file_path.unlink()
        logger.info(f"Deleted processed file: {file_path}")
```

---

## 2. Filtering Pipeline

### 2.1 MeSH Terms for Genetics Filtering

**Primary MeSH Descriptors (High Relevance):**

```python
# Primary MeSH terms - papers with these are almost certainly relevant
PRIMARY_MESH_TERMS = {
    # Genetics core
    "Genetics",
    "Genetic Variation",
    "Polymorphism, Single Nucleotide",
    "Polymorphism, Genetic",
    "Genotype",
    "Phenotype",
    "Alleles",
    "Gene Expression",
    "Gene Expression Regulation",

    # Specific genetic concepts
    "DNA Sequence Variation",
    "Genome-Wide Association Study",
    "Pharmacogenetics",
    "Pharmacogenomic Variants",
    "Genetic Predisposition to Disease",
    "Genetic Association Studies",

    # Molecular biology
    "DNA Methylation",
    "Epigenesis, Genetic",
    "MicroRNAs",
    "Transcription Factors",

    # Clinical genetics
    "Genetic Testing",
    "Genetic Counseling",
    "Precision Medicine",
    "Genetic Markers",
}

# Secondary MeSH terms - relevant when combined with genetics terms
SECONDARY_MESH_TERMS = {
    # Gene/protein related
    "Genes",
    "Proteins",
    "Enzymes",
    "Receptors",
    "Chromosomes",

    # Conditions often studied genetically
    "Neoplasms",
    "Cardiovascular Diseases",
    "Metabolic Diseases",
    "Autoimmune Diseases",
    "Neurodegenerative Diseases",

    # Pathways
    "Signal Transduction",
    "Metabolic Networks and Pathways",
    "Drug Interactions",
}

# MeSH subheadings that indicate genetic content
GENETIC_SUBHEADINGS = {
    "/genetics",
    "/metabolism",
    "/pharmacology",
    "/drug effects",
}
```

### 2.2 Keyword Filtering for SNPs, Genes, Pathways

```python
# scripts/etl/filter/genetics_filter.py

import re
from typing import Set, Optional
from dataclasses import dataclass

@dataclass
class FilterResult:
    is_relevant: bool
    confidence: float  # 0.0 - 1.0
    matched_terms: Set[str]
    matched_genes: Set[str]
    matched_snps: Set[str]
    filter_type: str  # 'mesh', 'keyword', 'entity'

class GeneticsFilter:
    """Filter papers for genetics relevance."""

    # SNP patterns
    SNP_PATTERNS = [
        re.compile(r'\brs\d{1,12}\b', re.IGNORECASE),           # rs1801133
        re.compile(r'\b[ACGT]>\b[ACGT]\b'),                      # A>G notation
        re.compile(r'\bp\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\b'),     # p.Ala222Val
        re.compile(r'\bc\.\d+[ACGT]>[ACGT]\b'),                  # c.677C>T
    ]

    # Gene symbol pattern (2-10 uppercase letters/numbers)
    GENE_PATTERN = re.compile(r'\b[A-Z][A-Z0-9]{1,9}\b')

    # Common non-gene abbreviations to exclude
    NON_GENE_TERMS = {
        'DNA', 'RNA', 'PCR', 'SNP', 'GWAS', 'BMI', 'LDL', 'HDL',
        'FDA', 'WHO', 'USA', 'UK', 'CI', 'OR', 'HR', 'RR', 'SD',
        'MRI', 'CT', 'PET', 'ECG', 'EEG', 'IV', 'IM', 'SC',
        'QTL', 'MAF', 'HWE', 'LD', 'PCA', 'FDR', 'ANOVA',
    }

    # Keywords indicating genetics content
    GENETICS_KEYWORDS = {
        # SNP/variant related
        'polymorphism', 'variant', 'mutation', 'allele', 'genotype',
        'haplotype', 'diplotype', 'heterozygous', 'homozygous',
        'minor allele', 'reference allele', 'wild type', 'wild-type',

        # Study types
        'genome-wide', 'gwas', 'association study', 'genetic association',
        'linkage', 'heritability', 'polygenic', 'mendelian randomization',

        # Pharmacogenomics
        'pharmacogenetic', 'pharmacogenomic', 'drug metabolism',
        'drug response', 'adverse drug reaction', 'dosing recommendation',
        'cyp2d6', 'cyp2c19', 'cyp3a4', 'vkorc1', 'dpyd', 'tpmt',

        # Molecular mechanisms
        'gene expression', 'transcription', 'translation',
        'methylation', 'epigenetic', 'splicing', 'promoter',

        # Pathways
        'pathway', 'signaling', 'metabolism', 'biosynthesis',
    }

    def __init__(self, known_genes: Set[str] = None):
        """
        Initialize filter with known gene symbols.

        Args:
            known_genes: Set of valid gene symbols (e.g., from NCBI Gene)
        """
        self.known_genes = known_genes or set()

    def filter_paper(self,
                     title: str,
                     abstract: str,
                     mesh_terms: list[str],
                     keywords: list[str] = None) -> FilterResult:
        """
        Determine if a paper is genetics-relevant.

        Returns FilterResult with confidence score and matched entities.
        """
        matched_terms = set()
        matched_genes = set()
        matched_snps = set()
        confidence = 0.0

        text = f"{title} {abstract}".lower()
        text_upper = f"{title} {abstract}"

        # 1. MeSH term matching (highest confidence)
        mesh_lower = {m.lower() for m in mesh_terms}

        for mesh in PRIMARY_MESH_TERMS:
            if mesh.lower() in mesh_lower:
                matched_terms.add(f"MeSH:{mesh}")
                confidence += 0.3

        for mesh in SECONDARY_MESH_TERMS:
            if mesh.lower() in mesh_lower:
                matched_terms.add(f"MeSH:{mesh}")
                confidence += 0.1

        # Check for genetic subheadings
        for mesh in mesh_terms:
            for subheading in GENETIC_SUBHEADINGS:
                if subheading in mesh.lower():
                    matched_terms.add(f"Subheading:{subheading}")
                    confidence += 0.15

        # 2. SNP extraction
        for pattern in self.SNP_PATTERNS:
            matches = pattern.findall(text_upper)
            for match in matches:
                matched_snps.add(match.lower())

        if matched_snps:
            confidence += min(0.3, len(matched_snps) * 0.1)

        # 3. Gene symbol extraction
        potential_genes = self.GENE_PATTERN.findall(text_upper)
        for gene in potential_genes:
            if gene not in self.NON_GENE_TERMS:
                if self.known_genes:
                    if gene in self.known_genes:
                        matched_genes.add(gene)
                else:
                    # Without gene list, be more conservative
                    if len(gene) >= 3 and len(gene) <= 6:
                        matched_genes.add(gene)

        if matched_genes:
            confidence += min(0.2, len(matched_genes) * 0.05)

        # 4. Keyword matching
        for keyword in self.GENETICS_KEYWORDS:
            if keyword in text:
                matched_terms.add(f"keyword:{keyword}")
                confidence += 0.1

        # Cap confidence at 1.0
        confidence = min(1.0, confidence)

        # Determine filter type
        if any('MeSH' in t for t in matched_terms):
            filter_type = 'mesh'
        elif matched_snps:
            filter_type = 'snp'
        elif matched_genes:
            filter_type = 'gene'
        else:
            filter_type = 'keyword'

        return FilterResult(
            is_relevant=confidence >= 0.3,  # Threshold for relevance
            confidence=confidence,
            matched_terms=matched_terms,
            matched_genes=matched_genes,
            matched_snps=matched_snps,
            filter_type=filter_type
        )

    def quick_filter(self, text: str) -> bool:
        """
        Fast preliminary filter before full processing.

        Use for initial pass to reduce processing load.
        """
        text_lower = text.lower()

        # Quick checks for high-confidence indicators
        quick_patterns = [
            r'\brs\d+\b',           # SNP rs numbers
            r'gwas\b',              # GWAS studies
            r'polymorphism',        # Polymorphism
            r'genetic variant',     # Genetic variant
            r'gene expression',     # Gene expression
            r'pharmacogene',        # Pharmacogenetics
        ]

        for pattern in quick_patterns:
            if re.search(pattern, text_lower):
                return True

        return False
```

### 2.3 Publication Type Filtering

```python
# Publication types to include
INCLUDE_PUBLICATION_TYPES = {
    "Journal Article",
    "Research Support, N.I.H., Extramural",
    "Research Support, Non-U.S. Gov't",
    "Research Support, U.S. Gov't, P.H.S.",
    "Multicenter Study",
    "Clinical Trial",
    "Randomized Controlled Trial",
    "Meta-Analysis",
    "Systematic Review",
    "Observational Study",
    "Comparative Study",
    "Validation Study",
    "Twin Study",
    "Genome-Wide Association Study",  # GWAS publication type
}

# Publication types to exclude (or deprioritize)
EXCLUDE_PUBLICATION_TYPES = {
    "Case Reports",         # Usually too specific
    "Letter",               # Low content
    "Editorial",            # Opinion
    "Comment",              # Response
    "News",                 # Not research
    "Biography",            # Not research
    "Historical Article",   # Limited relevance
    "Retracted Publication",  # Invalid
}

# Publication types to handle specially
SPECIAL_PUBLICATION_TYPES = {
    "Review": "include_if_genetics",      # Reviews are valuable for RAG
    "Practice Guideline": "high_priority", # Clinical guidelines
}
```

### 2.4 Estimated Yield

```
PubMed Total: ~36M papers
    │
    ├─► MeSH Filter (primary genetics terms)
    │   └─► ~8M papers (22%)
    │
    ├─► SNP/Gene Keyword Filter
    │   └─► ~5M papers (14%)
    │
    ├─► Combined Filter (MeSH OR Keywords)
    │   └─► ~10M papers (28%)
    │
    ├─► Confidence >= 0.3 (our threshold)
    │   └─► ~5M papers (14%)
    │
    ├─► Confidence >= 0.5 (high confidence)
    │   └─► ~3M papers (8%)
    │
    └─► Target: 1-3M highest relevance papers

Publication Year Distribution:
├── 2020-2024: ~1.5M papers (rapid growth in genetics)
├── 2015-2019: ~1.0M papers
├── 2010-2014: ~0.7M papers
└── Pre-2010:  ~0.8M papers (declining relevance)
```

---

## 3. Processing Pipeline

### 3.1 XML Parsing

```python
# scripts/etl/process/pubmed_parser.py

from lxml import etree
from typing import Generator, Optional
from dataclasses import dataclass, field
from pathlib import Path
import gzip

@dataclass
class PubMedArticle:
    """Parsed PubMed article structure."""
    pmid: int
    title: str
    abstract: str
    authors: list[dict] = field(default_factory=list)
    journal: str = ""
    publication_date: str = ""
    publication_types: list[str] = field(default_factory=list)
    mesh_terms: list[str] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    doi: Optional[str] = None
    pmc_id: Optional[str] = None
    language: str = "eng"

    # Extracted entities (filled later)
    extracted_genes: list[str] = field(default_factory=list)
    extracted_snps: list[str] = field(default_factory=list)

    # Relevance scoring (filled later)
    relevance_score: float = 0.0
    relevance_type: str = ""


class PubMedXMLParser:
    """
    Parse PubMed XML files using lxml for performance.

    lxml is ~10x faster than BeautifulSoup for XML parsing.
    Memory-efficient streaming parser for large files.
    """

    # XML namespaces (if needed)
    NAMESPACES = {}

    def parse_file(self, file_path: Path) -> Generator[PubMedArticle, None, None]:
        """
        Stream parse a PubMed XML file.

        Uses iterparse for memory efficiency - processes one article at a time
        without loading entire file into memory.
        """
        open_func = gzip.open if str(file_path).endswith('.gz') else open

        with open_func(file_path, 'rb') as f:
            # Use iterparse for memory-efficient streaming
            context = etree.iterparse(
                f,
                events=('end',),
                tag='PubmedArticle'
            )

            for event, elem in context:
                try:
                    article = self._parse_article(elem)
                    if article:
                        yield article
                except Exception as e:
                    pmid = elem.findtext('.//PMID')
                    print(f"Error parsing PMID {pmid}: {e}")
                finally:
                    # Clear element to free memory
                    elem.clear()
                    # Also clear parent references
                    while elem.getprevious() is not None:
                        del elem.getparent()[0]

    def _parse_article(self, elem: etree._Element) -> Optional[PubMedArticle]:
        """Parse a single PubmedArticle element."""

        # Get PMID
        pmid_elem = elem.find('.//PMID')
        if pmid_elem is None:
            return None
        pmid = int(pmid_elem.text)

        # Get article metadata
        article = elem.find('.//Article')
        if article is None:
            return None

        # Title
        title_elem = article.find('.//ArticleTitle')
        title = self._get_text_content(title_elem) if title_elem is not None else ""

        # Abstract
        abstract_elem = article.find('.//Abstract')
        abstract = self._parse_abstract(abstract_elem) if abstract_elem is not None else ""

        # Authors
        authors = self._parse_authors(article.find('.//AuthorList'))

        # Journal
        journal_elem = article.find('.//Journal/Title')
        journal = journal_elem.text if journal_elem is not None else ""

        # Publication date
        pub_date = self._parse_pub_date(article.find('.//ArticleDate') or
                                        article.find('.//Journal/JournalIssue/PubDate'))

        # Publication types
        pub_types = [
            pt.text for pt in article.findall('.//PublicationTypeList/PublicationType')
            if pt.text
        ]

        # MeSH terms
        mesh_list = elem.find('.//MeshHeadingList')
        mesh_terms = self._parse_mesh_terms(mesh_list)

        # Keywords
        keywords = [
            kw.text for kw in elem.findall('.//KeywordList/Keyword')
            if kw.text
        ]

        # DOI
        doi = None
        for aid in article.findall('.//ELocationID'):
            if aid.get('EIdType') == 'doi':
                doi = aid.text
                break

        # PMC ID
        pmc_id = None
        for aid in elem.findall('.//PubmedData/ArticleIdList/ArticleId'):
            if aid.get('IdType') == 'pmc':
                pmc_id = aid.text
                break

        # Language
        lang_elem = article.find('.//Language')
        language = lang_elem.text if lang_elem is not None else "eng"

        return PubMedArticle(
            pmid=pmid,
            title=title,
            abstract=abstract,
            authors=authors,
            journal=journal,
            publication_date=pub_date,
            publication_types=pub_types,
            mesh_terms=mesh_terms,
            keywords=keywords,
            doi=doi,
            pmc_id=pmc_id,
            language=language,
        )

    def _get_text_content(self, elem: etree._Element) -> str:
        """Get all text content including nested elements."""
        if elem is None:
            return ""
        return ''.join(elem.itertext()).strip()

    def _parse_abstract(self, abstract_elem: etree._Element) -> str:
        """Parse abstract, handling structured abstracts."""
        if abstract_elem is None:
            return ""

        parts = []
        for text in abstract_elem.findall('.//AbstractText'):
            label = text.get('Label', '')
            content = self._get_text_content(text)

            if label and content:
                parts.append(f"{label}: {content}")
            elif content:
                parts.append(content)

        return ' '.join(parts)

    def _parse_authors(self, author_list: Optional[etree._Element]) -> list[dict]:
        """Parse author information."""
        if author_list is None:
            return []

        authors = []
        for author in author_list.findall('.//Author'):
            lastname = author.findtext('LastName') or ''
            forename = author.findtext('ForeName') or author.findtext('Initials') or ''

            if lastname or forename:
                authors.append({
                    'lastname': lastname,
                    'forename': forename,
                    'full_name': f"{forename} {lastname}".strip()
                })

        return authors

    def _parse_pub_date(self, date_elem: Optional[etree._Element]) -> str:
        """Parse publication date into ISO format."""
        if date_elem is None:
            return ""

        year = date_elem.findtext('Year') or ''
        month = date_elem.findtext('Month') or '01'
        day = date_elem.findtext('Day') or '01'

        # Handle text months
        month_map = {
            'Jan': '01', 'Feb': '02', 'Mar': '03', 'Apr': '04',
            'May': '05', 'Jun': '06', 'Jul': '07', 'Aug': '08',
            'Sep': '09', 'Oct': '10', 'Nov': '11', 'Dec': '12'
        }

        if month in month_map:
            month = month_map[month]
        elif not month.isdigit():
            month = '01'

        month = month.zfill(2)
        day = day.zfill(2)

        if year:
            return f"{year}-{month}-{day}"
        return ""

    def _parse_mesh_terms(self, mesh_list: Optional[etree._Element]) -> list[str]:
        """Parse MeSH headings with subheadings."""
        if mesh_list is None:
            return []

        terms = []
        for heading in mesh_list.findall('.//MeshHeading'):
            descriptor = heading.find('DescriptorName')
            if descriptor is None:
                continue

            main_term = descriptor.text
            terms.append(main_term)

            # Add qualifiers/subheadings
            for qualifier in heading.findall('.//QualifierName'):
                if qualifier.text:
                    terms.append(f"{main_term}/{qualifier.text}")

        return terms
```

### 3.2 Text Cleaning and Normalization

```python
# scripts/etl/process/text_cleaner.py

import re
import unicodedata
from typing import Optional

class TextCleaner:
    """Clean and normalize text for embedding generation."""

    # Patterns to clean
    CITATION_PATTERN = re.compile(r'\[\d+(?:,\s*\d+)*\]|\(\d+(?:,\s*\d+)*\)')
    URL_PATTERN = re.compile(r'https?://\S+|www\.\S+')
    EMAIL_PATTERN = re.compile(r'\S+@\S+\.\S+')
    MULTIPLE_SPACES = re.compile(r'\s+')
    SPECIAL_CHARS = re.compile(r'[^\w\s\-\.,;:!?\(\)\'\"]+')

    @classmethod
    def clean_for_embedding(cls, text: str) -> str:
        """
        Clean text for embedding generation.

        - Removes citations [1,2,3]
        - Normalizes unicode
        - Removes URLs/emails
        - Normalizes whitespace
        - Preserves scientific notation
        """
        if not text:
            return ""

        # Unicode normalization
        text = unicodedata.normalize('NFKD', text)

        # Remove citations
        text = cls.CITATION_PATTERN.sub('', text)

        # Remove URLs and emails
        text = cls.URL_PATTERN.sub('', text)
        text = cls.EMAIL_PATTERN.sub('', text)

        # Normalize whitespace
        text = cls.MULTIPLE_SPACES.sub(' ', text)

        return text.strip()

    @classmethod
    def clean_for_display(cls, text: str) -> str:
        """
        Clean text for display while preserving formatting.
        """
        if not text:
            return ""

        # Just normalize unicode and whitespace
        text = unicodedata.normalize('NFKD', text)
        text = cls.MULTIPLE_SPACES.sub(' ', text)

        return text.strip()

    @classmethod
    def truncate_for_embedding(cls, text: str, max_tokens: int = 256) -> str:
        """
        Truncate text to fit embedding model's context window.

        all-MiniLM-L6-v2 has 256 token limit.
        Rough estimate: 1 token ~= 4 characters for English.
        """
        max_chars = max_tokens * 4  # Conservative estimate

        if len(text) <= max_chars:
            return text

        # Try to truncate at sentence boundary
        truncated = text[:max_chars]
        last_period = truncated.rfind('. ')

        if last_period > max_chars * 0.7:  # If we can keep at least 70%
            return truncated[:last_period + 1]

        return truncated.strip()
```

### 3.3 Entity Extraction

```python
# scripts/etl/process/entity_extractor.py

import re
from typing import Set, Tuple
from dataclasses import dataclass

@dataclass
class ExtractedEntities:
    """Container for extracted biomedical entities."""
    genes: Set[str]
    snps: Set[str]
    drugs: Set[str]
    diseases: Set[str]
    pathways: Set[str]


class BioEntityExtractor:
    """
    Extract biomedical entities from text.

    Uses pattern-based extraction for speed.
    For production, consider NER models like BioBERT or PubMedBERT.
    """

    # SNP patterns
    SNP_RS = re.compile(r'\b(rs\d{1,12})\b', re.IGNORECASE)
    SNP_HGVS_PROTEIN = re.compile(r'\b(p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2})\b')
    SNP_HGVS_CODING = re.compile(r'\b(c\.\d+[ACGT]>[ACGT])\b')

    # Gene patterns (validated against known genes)
    GENE_PATTERN = re.compile(r'\b([A-Z][A-Z0-9]{1,9})\b')

    # Drug patterns (common suffixes)
    DRUG_SUFFIXES = ('mab', 'nib', 'zole', 'pril', 'sartan', 'statin',
                     'vir', 'cycline', 'mycin', 'cillin', 'pam', 'lol')

    def __init__(self,
                 known_genes: Set[str] = None,
                 known_drugs: Set[str] = None):
        self.known_genes = known_genes or set()
        self.known_drugs = known_drugs or set()

        # Common abbreviations that look like genes but aren't
        self.non_genes = {
            'DNA', 'RNA', 'PCR', 'SNP', 'GWAS', 'BMI', 'LDL', 'HDL',
            'FDA', 'WHO', 'USA', 'UK', 'EU', 'CI', 'OR', 'HR', 'RR',
            'SD', 'SE', 'NS', 'NA', 'ND', 'NR', 'vs', 'mg', 'kg',
            'mL', 'dL', 'mmol', 'nmol', 'pmol', 'ng', 'pg', 'IU',
        }

    def extract_all(self, text: str) -> ExtractedEntities:
        """Extract all entity types from text."""
        return ExtractedEntities(
            genes=self.extract_genes(text),
            snps=self.extract_snps(text),
            drugs=self.extract_drugs(text),
            diseases=set(),  # Requires NER for accuracy
            pathways=set(),  # Requires knowledge base lookup
        )

    def extract_snps(self, text: str) -> Set[str]:
        """
        Extract SNP identifiers.

        Extracts:
        - rs numbers (rs1801133)
        - HGVS protein notation (p.Ala222Val)
        - HGVS coding notation (c.677C>T)
        """
        snps = set()

        # rs numbers (normalized to lowercase)
        for match in self.SNP_RS.finditer(text):
            snps.add(match.group(1).lower())

        # HGVS protein notation
        for match in self.SNP_HGVS_PROTEIN.finditer(text):
            snps.add(match.group(1))

        # HGVS coding notation
        for match in self.SNP_HGVS_CODING.finditer(text):
            snps.add(match.group(1))

        return snps

    def extract_genes(self, text: str) -> Set[str]:
        """
        Extract gene symbols.

        If known_genes provided, validates against it.
        Otherwise uses heuristics.
        """
        genes = set()

        for match in self.GENE_PATTERN.finditer(text):
            symbol = match.group(1)

            # Skip common non-genes
            if symbol in self.non_genes:
                continue

            # Validate against known genes if available
            if self.known_genes:
                if symbol in self.known_genes:
                    genes.add(symbol)
            else:
                # Heuristics for likely gene symbols
                # - 2-6 characters
                # - Contains at least one letter
                # - Common gene naming patterns
                if 2 <= len(symbol) <= 6:
                    genes.add(symbol)

        return genes

    def extract_drugs(self, text: str) -> Set[str]:
        """
        Extract drug names.

        Uses suffix patterns and known drug list.
        """
        drugs = set()

        # Check known drugs
        text_lower = text.lower()
        for drug in self.known_drugs:
            if drug.lower() in text_lower:
                drugs.add(drug)

        # Pattern-based extraction
        words = re.findall(r'\b\w+\b', text)
        for word in words:
            word_lower = word.lower()
            for suffix in self.DRUG_SUFFIXES:
                if word_lower.endswith(suffix) and len(word) > len(suffix) + 2:
                    drugs.add(word)
                    break

        return drugs


class SectionIdentifier:
    """
    Identify sections in structured abstracts.

    Useful for prioritizing certain sections for embedding.
    """

    # Standard IMRAD sections
    SECTION_PATTERNS = {
        'background': ['background', 'introduction', 'context', 'rationale'],
        'objective': ['objective', 'objectives', 'aim', 'aims', 'purpose'],
        'methods': ['method', 'methods', 'design', 'materials', 'patients',
                    'participants', 'setting', 'intervention'],
        'results': ['result', 'results', 'findings', 'outcome', 'outcomes'],
        'conclusion': ['conclusion', 'conclusions', 'interpretation',
                       'implications', 'summary', 'significance'],
    }

    @classmethod
    def identify_sections(cls, abstract: str) -> dict:
        """
        Parse structured abstract into sections.

        Returns dict mapping section type to content.
        """
        sections = {}

        # Try to split on common patterns
        # Pattern: "BACKGROUND: text METHODS: text RESULTS: text"
        pattern = r'([A-Z][A-Z\s]+):\s*'
        parts = re.split(pattern, abstract)

        if len(parts) > 1:
            # Parts alternate: [before_first_label, label1, content1, label2, content2, ...]
            for i in range(1, len(parts) - 1, 2):
                label = parts[i].strip().lower()
                content = parts[i + 1].strip()

                # Map to standard section types
                for section_type, keywords in cls.SECTION_PATTERNS.items():
                    if any(kw in label for kw in keywords):
                        sections[section_type] = content
                        break

        return sections

    @classmethod
    def get_key_sections(cls, abstract: str) -> str:
        """
        Extract key sections for embedding.

        Prioritizes: Results > Conclusion > Methods
        """
        sections = cls.identify_sections(abstract)

        # Priority order for embedding
        priority = ['results', 'conclusion', 'objective', 'background']

        key_content = []
        for section in priority:
            if section in sections:
                key_content.append(sections[section])

        if key_content:
            return ' '.join(key_content)

        return abstract  # Return full abstract if no sections found
```

---

## 4. Embedding Generation

### 4.1 Model Selection

| Model | Dimensions | Speed (CPU) | Speed (GPU) | Quality | Size |
|-------|------------|-------------|-------------|---------|------|
| **all-MiniLM-L6-v2** | 384 | ~1,000/sec | ~15,000/sec | Excellent | 80 MB |
| all-mpnet-base-v2 | 768 | ~400/sec | ~8,000/sec | Better | 420 MB |
| PubMedBERT | 768 | ~300/sec | ~5,000/sec | Biomedical | 440 MB |
| BioLinkBERT | 768 | ~300/sec | ~5,000/sec | Entity linking | 440 MB |

**Recommendation: all-MiniLM-L6-v2**

Reasons:
1. **RuVector optimized for 384d** - Built-in tiered compression
2. **4x faster** than 768d alternatives
3. **Quality sufficient** - 0.82 vs 0.84 on STS benchmarks
4. **50% smaller storage** - Critical for cost-sensitive deployment
5. **RuVector ONNX runtime** - No external API needed, runs locally

### 4.2 Chunking Strategy

```python
# scripts/etl/embeddings/chunker.py

from dataclasses import dataclass
from typing import List, Optional
import re

@dataclass
class TextChunk:
    """A chunk of text for embedding."""
    text: str
    chunk_index: int
    total_chunks: int
    source_field: str  # 'title', 'abstract', 'full_text'
    metadata: dict


class ArticleChunker:
    """
    Chunk articles for embedding generation.

    Strategy for research papers:
    1. Title alone (always)
    2. Title + Abstract combined (if fits)
    3. Abstract sections separately (if structured)
    4. Full text paragraphs (if available)
    """

    # all-MiniLM-L6-v2 context window
    MAX_TOKENS = 256
    MAX_CHARS = 1000  # Conservative: 256 * 4

    # Overlap for context preservation
    OVERLAP_CHARS = 100

    def __init__(self, max_chars: int = MAX_CHARS):
        self.max_chars = max_chars

    def chunk_article(self,
                      pmid: int,
                      title: str,
                      abstract: str,
                      full_text: Optional[str] = None) -> List[TextChunk]:
        """
        Generate chunks for a single article.

        For most papers (abstract only), generates 1-3 chunks.
        """
        chunks = []

        # Chunk 1: Title alone (for title-based search)
        chunks.append(TextChunk(
            text=title,
            chunk_index=0,
            total_chunks=0,  # Will update
            source_field='title',
            metadata={'pmid': pmid}
        ))

        # Chunk 2: Title + Abstract (primary search chunk)
        title_abstract = f"{title}. {abstract}"
        if len(title_abstract) <= self.max_chars:
            chunks.append(TextChunk(
                text=title_abstract,
                chunk_index=1,
                total_chunks=0,
                source_field='title_abstract',
                metadata={'pmid': pmid}
            ))
        else:
            # Abstract is too long, need to chunk it
            abstract_chunks = self._chunk_text(
                abstract,
                prefix=title[:100] + ". "  # Include title prefix for context
            )
            for i, chunk_text in enumerate(abstract_chunks):
                chunks.append(TextChunk(
                    text=chunk_text,
                    chunk_index=len(chunks),
                    total_chunks=0,
                    source_field='abstract',
                    metadata={'pmid': pmid, 'abstract_chunk': i}
                ))

        # Chunk 3+: Full text (if available)
        if full_text:
            full_text_chunks = self._chunk_text(full_text)
            for i, chunk_text in enumerate(full_text_chunks):
                chunks.append(TextChunk(
                    text=chunk_text,
                    chunk_index=len(chunks),
                    total_chunks=0,
                    source_field='full_text',
                    metadata={'pmid': pmid, 'full_text_chunk': i}
                ))

        # Update total_chunks
        for chunk in chunks:
            chunk.total_chunks = len(chunks)

        return chunks

    def _chunk_text(self, text: str, prefix: str = "") -> List[str]:
        """
        Chunk text with overlap.

        Tries to break at sentence boundaries.
        """
        chunks = []
        prefix_len = len(prefix)
        available_chars = self.max_chars - prefix_len

        # Split into sentences
        sentences = re.split(r'(?<=[.!?])\s+', text)

        current_chunk = prefix
        for sentence in sentences:
            if len(current_chunk) + len(sentence) + 1 <= self.max_chars:
                current_chunk += sentence + " "
            else:
                if current_chunk.strip():
                    chunks.append(current_chunk.strip())

                # Start new chunk with overlap
                if chunks:
                    # Get last 100 chars for overlap
                    overlap = current_chunk[-self.OVERLAP_CHARS:] if len(current_chunk) > self.OVERLAP_CHARS else ""
                    current_chunk = overlap + sentence + " "
                else:
                    current_chunk = prefix + sentence + " "

        if current_chunk.strip():
            chunks.append(current_chunk.strip())

        return chunks

    def chunk_for_rag(self,
                      pmid: int,
                      title: str,
                      abstract: str) -> TextChunk:
        """
        Generate single optimal chunk for RAG retrieval.

        Combines title + key abstract sections within token limit.
        """
        # For most abstracts, title + abstract fits
        combined = f"{title}. {abstract}"

        if len(combined) <= self.max_chars:
            return TextChunk(
                text=combined,
                chunk_index=0,
                total_chunks=1,
                source_field='combined',
                metadata={'pmid': pmid}
            )

        # If too long, prioritize results/conclusion
        from .entity_extractor import SectionIdentifier
        key_text = SectionIdentifier.get_key_sections(abstract)

        combined = f"{title}. {key_text}"

        if len(combined) <= self.max_chars:
            return TextChunk(
                text=combined,
                chunk_index=0,
                total_chunks=1,
                source_field='key_sections',
                metadata={'pmid': pmid}
            )

        # Last resort: truncate
        return TextChunk(
            text=combined[:self.max_chars],
            chunk_index=0,
            total_chunks=1,
            source_field='truncated',
            metadata={'pmid': pmid}
        )
```

### 4.3 Batch Embedding Generation

```python
# scripts/etl/embeddings/generator.py

from pathlib import Path
from typing import List, Generator, Optional
import numpy as np
import json
from dataclasses import asdict

# RuVector's built-in ONNX embedding
# Alternative: sentence-transformers for more model options

class EmbeddingGenerator:
    """
    Generate embeddings using RuVector's ONNX runtime.

    RuVector includes all-MiniLM-L6-v2 by default.
    No external API calls needed.
    """

    def __init__(self,
                 model_name: str = "all-MiniLM-L6-v2",
                 batch_size: int = 64,
                 device: str = "cpu"):
        """
        Initialize embedding generator.

        Args:
            model_name: Model to use (RuVector supports MiniLM by default)
            batch_size: Batch size for embedding generation
            device: 'cpu' or 'cuda' (if GPU available)
        """
        self.model_name = model_name
        self.batch_size = batch_size
        self.device = device
        self.dimension = 384

        # For production: use RuVector's built-in ONNX
        # For development: can use sentence-transformers
        try:
            from sentence_transformers import SentenceTransformer
            self.model = SentenceTransformer(model_name, device=device)
            self.use_sentence_transformers = True
        except ImportError:
            # Fall back to RuVector ONNX (Node.js)
            self.model = None
            self.use_sentence_transformers = False

    def generate_batch(self, texts: List[str]) -> np.ndarray:
        """
        Generate embeddings for a batch of texts.

        Returns numpy array of shape (len(texts), 384).
        """
        if self.use_sentence_transformers:
            embeddings = self.model.encode(
                texts,
                batch_size=self.batch_size,
                show_progress_bar=False,
                convert_to_numpy=True,
                normalize_embeddings=True
            )
            return embeddings
        else:
            raise RuntimeError("Please install sentence-transformers or use RuVector Node.js")

    def process_articles(self,
                         articles: Generator,
                         output_path: Path,
                         chunker) -> int:
        """
        Process articles and save embeddings.

        Streams processing to handle large datasets.
        """
        batch_texts = []
        batch_metadata = []
        total_processed = 0

        with open(output_path, 'w') as f:
            for article in articles:
                # Generate chunks for this article
                chunk = chunker.chunk_for_rag(
                    pmid=article.pmid,
                    title=article.title,
                    abstract=article.abstract
                )

                batch_texts.append(chunk.text)
                batch_metadata.append({
                    'pmid': article.pmid,
                    'title': article.title[:200],  # Truncate for storage
                    'chunk_index': chunk.chunk_index,
                    'source_field': chunk.source_field,
                })

                # Process batch when full
                if len(batch_texts) >= self.batch_size:
                    embeddings = self.generate_batch(batch_texts)

                    for i, (emb, meta) in enumerate(zip(embeddings, batch_metadata)):
                        record = {
                            **meta,
                            'embedding': emb.tolist()
                        }
                        f.write(json.dumps(record) + '\n')

                    total_processed += len(batch_texts)
                    print(f"Processed {total_processed} articles...")

                    batch_texts = []
                    batch_metadata = []

            # Process remaining batch
            if batch_texts:
                embeddings = self.generate_batch(batch_texts)
                for emb, meta in zip(embeddings, batch_metadata):
                    record = {**meta, 'embedding': emb.tolist()}
                    f.write(json.dumps(record) + '\n')
                total_processed += len(batch_texts)

        return total_processed


def estimate_processing_time(num_articles: int, device: str = 'cpu') -> dict:
    """
    Estimate processing time for embedding generation.

    Returns dict with time estimates.
    """
    if device == 'cpu':
        articles_per_second = 1000
    else:  # GPU
        articles_per_second = 15000

    total_seconds = num_articles / articles_per_second

    return {
        'articles': num_articles,
        'device': device,
        'articles_per_second': articles_per_second,
        'total_seconds': total_seconds,
        'total_minutes': total_seconds / 60,
        'total_hours': total_seconds / 3600,
    }

# Example estimates:
# 1M articles, CPU: ~17 minutes
# 1M articles, GPU: ~1 minute
# 5M articles, CPU: ~1.4 hours
# 5M articles, GPU: ~5.5 minutes
```

### 4.4 GPU vs CPU Considerations

| Consideration | CPU | GPU |
|---------------|-----|-----|
| **Speed** | ~1,000/sec | ~15,000/sec |
| **Cost (cloud)** | ~$0.05/hour | ~$0.50/hour |
| **Memory** | 4 GB sufficient | 8 GB VRAM recommended |
| **Setup** | Simple | Requires CUDA |
| **1M articles** | 17 minutes | 1 minute |
| **10M articles** | 2.8 hours | 11 minutes |

**Recommendation for bootstrapped startup:**

1. **Initial load**: Rent GPU instance for 1-2 hours (~$1-2)
2. **Daily updates**: CPU is fine (~1,000 new papers = 1 second)
3. **RuVector ONNX**: Use built-in CPU runtime (no setup)

---

## 5. Storage Pipeline

### 5.1 RuVector Ingestion

```python
# scripts/etl/load/ruvector_loader.py

import json
from pathlib import Path
from typing import List, Dict, Any
import subprocess

class RuVectorLoader:
    """
    Load embeddings into RuVector.

    Uses RuVector's CLI for bulk import.
    Alternative: Use RuVector Node.js SDK directly.
    """

    def __init__(self, collection: str = "articles"):
        self.collection = collection

    def prepare_import_file(self,
                            embeddings_path: Path,
                            output_path: Path,
                            batch_size: int = 10000) -> List[Path]:
        """
        Prepare embeddings for RuVector import.

        Converts JSONL to RuVector's expected format.
        Splits into batches for memory efficiency.
        """
        batch_files = []
        current_batch = []
        batch_num = 0

        with open(embeddings_path, 'r') as f:
            for line in f:
                record = json.loads(line)

                # Convert to RuVector format
                ruvector_record = {
                    'id': f"pmid:{record['pmid']}",
                    'embedding': record['embedding'],
                    'metadata': {
                        'pmid': record['pmid'],
                        'title': record['title'],
                        'source_field': record['source_field'],
                    }
                }

                current_batch.append(ruvector_record)

                if len(current_batch) >= batch_size:
                    batch_file = output_path / f"batch_{batch_num:04d}.json"
                    with open(batch_file, 'w') as bf:
                        json.dump(current_batch, bf)
                    batch_files.append(batch_file)
                    current_batch = []
                    batch_num += 1

        # Write remaining
        if current_batch:
            batch_file = output_path / f"batch_{batch_num:04d}.json"
            with open(batch_file, 'w') as bf:
                json.dump(current_batch, bf)
            batch_files.append(batch_file)

        return batch_files

    def import_batch(self, batch_file: Path) -> bool:
        """
        Import a batch file using RuVector CLI.
        """
        try:
            result = subprocess.run([
                'npx', 'ruvector', 'import',
                '--collection', self.collection,
                '--file', str(batch_file),
            ], capture_output=True, text=True)

            if result.returncode != 0:
                print(f"Import error: {result.stderr}")
                return False

            return True
        except Exception as e:
            print(f"Import failed: {e}")
            return False

    def create_graph_relationships(self, article_entities: List[Dict]):
        """
        Create graph relationships in RuVector.

        Links articles to genes, SNPs, pathways.
        """
        for article in article_entities:
            pmid = article['pmid']

            # Create CITED_IN relationships
            for gene in article.get('genes', []):
                self._create_relationship(
                    f"gene:{gene}",
                    f"pmid:{pmid}",
                    "CITED_IN"
                )

            for snp in article.get('snps', []):
                self._create_relationship(
                    f"snp:{snp}",
                    f"pmid:{pmid}",
                    "CITED_IN"
                )

    def _create_relationship(self, from_id: str, to_id: str, rel_type: str):
        """Create a single relationship in RuVector."""
        # Implementation depends on RuVector SDK
        pass
```

### 5.2 Supabase Storage

```python
# scripts/etl/load/supabase_loader.py

from supabase import create_client, Client
from typing import List, Dict, Any
import os

class SupabaseLoader:
    """
    Load article metadata into Supabase PostgreSQL.

    Vector embeddings go to RuVector.
    Metadata and full records go to Supabase.
    """

    def __init__(self):
        url = os.environ['SUPABASE_URL']
        key = os.environ['SUPABASE_KEY']
        self.client: Client = create_client(url, key)

    def upsert_articles(self, articles: List[Dict[str, Any]]):
        """
        Upsert article metadata to Supabase.

        Schema:
        - pmid (PK)
        - title
        - abstract (truncated for storage)
        - authors (JSONB)
        - journal
        - publication_date
        - doi
        - pmc_id
        - mesh_terms (array)
        - relevance_score
        - extracted_genes (array)
        - extracted_snps (array)
        - created_at
        - updated_at
        """
        # Prepare records
        records = []
        for article in articles:
            records.append({
                'pmid': article['pmid'],
                'title': article['title'][:500],
                'abstract': article['abstract'][:2000],  # Truncate for storage
                'authors': article.get('authors', []),
                'journal': article.get('journal', ''),
                'publication_date': article.get('publication_date'),
                'doi': article.get('doi'),
                'pmc_id': article.get('pmc_id'),
                'mesh_terms': article.get('mesh_terms', []),
                'relevance_score': article.get('relevance_score', 0),
                'extracted_genes': list(article.get('extracted_genes', [])),
                'extracted_snps': list(article.get('extracted_snps', [])),
            })

        # Batch upsert
        batch_size = 1000
        for i in range(0, len(records), batch_size):
            batch = records[i:i + batch_size]
            self.client.table('research_articles').upsert(
                batch,
                on_conflict='pmid'
            ).execute()
            print(f"Upserted {min(i + batch_size, len(records))}/{len(records)} articles")

    def create_article_entity_links(self, pmid: int, genes: List[str], snps: List[str]):
        """
        Create links between articles and extracted entities.

        Tables:
        - article_gene_links (pmid, gene_symbol)
        - article_snp_links (pmid, rs_number)
        """
        # Gene links
        for gene in genes:
            self.client.table('article_gene_links').upsert({
                'pmid': pmid,
                'gene_symbol': gene,
            }, on_conflict='pmid,gene_symbol').execute()

        # SNP links
        for snp in snps:
            self.client.table('article_snp_links').upsert({
                'pmid': pmid,
                'rs_number': snp,
            }, on_conflict='pmid,rs_number').execute()
```

### 5.3 Deduplication Strategy

```python
# scripts/etl/dedup/deduplicator.py

import hashlib
from typing import Set, Optional
from dataclasses import dataclass

@dataclass
class DedupResult:
    is_duplicate: bool
    duplicate_of: Optional[int]  # PMID of original
    reason: str


class ArticleDeduplicator:
    """
    Detect and handle duplicate articles.

    PubMed has some duplicates from:
    - Different publication stages (ahead of print vs final)
    - Corrections and errata
    - Republications
    """

    def __init__(self):
        # Cache of seen titles (normalized)
        self.title_hashes: dict[str, int] = {}

        # Cache of seen DOIs
        self.seen_dois: dict[str, int] = {}

        # Cache of seen content hashes
        self.content_hashes: dict[str, int] = {}

    def check_duplicate(self,
                        pmid: int,
                        title: str,
                        abstract: str,
                        doi: Optional[str] = None) -> DedupResult:
        """
        Check if article is a duplicate.

        Returns DedupResult with duplicate status and reason.
        """
        # 1. Check DOI (most reliable)
        if doi:
            if doi in self.seen_dois:
                return DedupResult(
                    is_duplicate=True,
                    duplicate_of=self.seen_dois[doi],
                    reason="duplicate_doi"
                )
            self.seen_dois[doi] = pmid

        # 2. Check title hash (for same article different formats)
        title_normalized = self._normalize_title(title)
        title_hash = hashlib.md5(title_normalized.encode()).hexdigest()

        if title_hash in self.title_hashes:
            return DedupResult(
                is_duplicate=True,
                duplicate_of=self.title_hashes[title_hash],
                reason="duplicate_title"
            )
        self.title_hashes[title_hash] = pmid

        # 3. Check content hash (for near-duplicates)
        content = f"{title_normalized} {self._normalize_text(abstract[:500])}"
        content_hash = hashlib.md5(content.encode()).hexdigest()

        if content_hash in self.content_hashes:
            return DedupResult(
                is_duplicate=True,
                duplicate_of=self.content_hashes[content_hash],
                reason="duplicate_content"
            )
        self.content_hashes[content_hash] = pmid

        return DedupResult(
            is_duplicate=False,
            duplicate_of=None,
            reason="unique"
        )

    def _normalize_title(self, title: str) -> str:
        """Normalize title for comparison."""
        import re
        # Lowercase, remove punctuation, normalize whitespace
        title = title.lower()
        title = re.sub(r'[^\w\s]', '', title)
        title = ' '.join(title.split())
        return title

    def _normalize_text(self, text: str) -> str:
        """Normalize text for hashing."""
        import re
        text = text.lower()
        text = re.sub(r'[^\w\s]', '', text)
        text = ' '.join(text.split())
        return text
```

---

## 6. Update Pipeline

### 6.1 Daily/Weekly Update Strategy

```python
# scripts/etl/updates/update_manager.py

from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
import json

class UpdateManager:
    """
    Manage incremental updates to the paper database.

    Strategy:
    - Daily: Download and process new papers
    - Weekly: Full quality check and reindexing
    - Monthly: Full refresh of high-value papers
    """

    def __init__(self, state_dir: Path):
        self.state_dir = state_dir
        self.state_file = state_dir / "update_state.json"
        self.state = self._load_state()

    def _load_state(self) -> dict:
        if self.state_file.exists():
            with open(self.state_file) as f:
                return json.load(f)
        return {
            "last_daily_update": None,
            "last_weekly_check": None,
            "last_monthly_refresh": None,
            "last_processed_update_file": None,
            "total_papers": 0,
            "update_history": []
        }

    def _save_state(self):
        with open(self.state_file, 'w') as f:
            json.dump(self.state, f, indent=2, default=str)

    def run_daily_update(self):
        """
        Process new papers from daily update files.

        Typical volume: 1,000-2,000 new papers/day
        Processing time: ~1-2 minutes
        """
        from ..extract.pubmed_downloader import PubMedDownloader
        from ..process.pubmed_parser import PubMedXMLParser
        from ..filter.genetics_filter import GeneticsFilter

        downloader = PubMedDownloader(self.state_dir / "downloads")
        parser = PubMedXMLParser()
        filter_ = GeneticsFilter()

        # Download updates since last run
        last_update = self.state.get("last_daily_update")
        if last_update:
            since_date = datetime.fromisoformat(last_update)
        else:
            since_date = datetime.now() - timedelta(days=7)

        new_papers = 0
        relevant_papers = 0

        for update_file in downloader.download_updates(since_date):
            for article in parser.parse_file(update_file):
                new_papers += 1

                # Filter for genetics relevance
                result = filter_.filter_paper(
                    title=article.title,
                    abstract=article.abstract,
                    mesh_terms=article.mesh_terms
                )

                if result.is_relevant:
                    relevant_papers += 1
                    # Process and store
                    self._process_new_article(article, result)

        # Update state
        self.state["last_daily_update"] = datetime.now().isoformat()
        self.state["total_papers"] += relevant_papers
        self.state["update_history"].append({
            "date": datetime.now().isoformat(),
            "type": "daily",
            "new_papers": new_papers,
            "relevant_papers": relevant_papers,
        })
        self._save_state()

        return {"new": new_papers, "relevant": relevant_papers}

    def run_weekly_check(self):
        """
        Weekly quality check and maintenance.

        - Check for retractions
        - Check for corrections
        - Re-index if needed
        - Clean up old files
        """
        retracted = self._check_retractions()
        corrected = self._check_corrections()

        self.state["last_weekly_check"] = datetime.now().isoformat()
        self._save_state()

        return {"retracted": retracted, "corrected": corrected}

    def _process_new_article(self, article, filter_result):
        """Process a single new article."""
        # Extract entities
        from ..process.entity_extractor import BioEntityExtractor
        extractor = BioEntityExtractor()
        entities = extractor.extract_all(f"{article.title} {article.abstract}")

        article.extracted_genes = list(entities.genes)
        article.extracted_snps = list(entities.snps)
        article.relevance_score = filter_result.confidence

        # Generate embedding
        # Store in RuVector and Supabase
        # (Implementation depends on batch vs streaming approach)
        pass

    def _check_retractions(self) -> int:
        """Check for retracted papers and mark them."""
        # Query PubMed for retractions
        # Update our database to mark/remove retracted papers
        return 0

    def _check_corrections(self) -> int:
        """Check for corrections and update records."""
        return 0
```

### 6.2 Handling Retractions and Corrections

```python
# scripts/etl/updates/retraction_handler.py

from typing import List, Set
import aiohttp
import asyncio

class RetractionHandler:
    """
    Handle retracted and corrected papers.

    PubMed marks retractions with:
    - Publication type: "Retracted Publication"
    - Separate retraction notices
    """

    RETRACTION_API = "https://api.retractiondatabase.org/api/v1/papers"

    async def get_retracted_pmids(self, since_date: str) -> Set[int]:
        """
        Get recently retracted PMIDs.

        Uses Retraction Watch database API.
        """
        async with aiohttp.ClientSession() as session:
            async with session.get(
                self.RETRACTION_API,
                params={"since": since_date, "format": "json"}
            ) as response:
                data = await response.json()
                return {int(p['pmid']) for p in data.get('papers', []) if p.get('pmid')}

    def mark_retracted(self, pmids: List[int]):
        """
        Mark papers as retracted in our database.

        Options:
        1. Delete entirely (data integrity)
        2. Mark as retracted but keep (historical record)
        3. Move to separate collection (searchable but flagged)

        Recommendation: Option 2 (mark but keep)
        """
        # Supabase update
        from supabase import create_client
        import os

        client = create_client(os.environ['SUPABASE_URL'], os.environ['SUPABASE_KEY'])

        for pmid in pmids:
            client.table('research_articles').update({
                'is_retracted': True,
                'retracted_at': 'now()'
            }).eq('pmid', pmid).execute()

        # RuVector: Add retraction flag to metadata
        # (Vectors kept for similarity, but flagged in results)


class CorrectionHandler:
    """
    Handle corrected papers.

    PubMed issues corrections as:
    - "Erratum" publication type
    - "Correction" publication type
    - Linked to original via related articles
    """

    def apply_corrections(self, corrections: List[dict]):
        """
        Apply corrections to existing papers.

        Strategy:
        1. Keep original record
        2. Add correction note
        3. Optionally re-embed if content significantly changed
        """
        for correction in corrections:
            original_pmid = correction['original_pmid']
            correction_text = correction['correction_text']

            # Update Supabase
            self._update_with_correction(original_pmid, correction_text)

            # Re-embed if abstract changed significantly
            if correction.get('abstract_changed'):
                self._regenerate_embedding(original_pmid)
```

### 6.3 Version Tracking

```sql
-- Supabase schema for version tracking

-- Article versions table
CREATE TABLE article_versions (
    id SERIAL PRIMARY KEY,
    pmid INTEGER NOT NULL,
    version INTEGER NOT NULL DEFAULT 1,

    -- Content at this version
    title TEXT,
    abstract TEXT,

    -- Change metadata
    change_type VARCHAR(50),  -- 'initial', 'correction', 'update'
    change_description TEXT,

    -- Timestamps
    created_at TIMESTAMP DEFAULT NOW(),

    UNIQUE(pmid, version)
);

-- Retraction tracking
CREATE TABLE retractions (
    id SERIAL PRIMARY KEY,
    pmid INTEGER NOT NULL UNIQUE,

    retraction_date DATE,
    retraction_notice_pmid INTEGER,
    reason TEXT,

    -- Our response
    action_taken VARCHAR(50),  -- 'flagged', 'removed', 'kept'
    action_date TIMESTAMP DEFAULT NOW()
);

-- Update log
CREATE TABLE etl_update_log (
    id SERIAL PRIMARY KEY,
    update_type VARCHAR(50),  -- 'daily', 'weekly', 'monthly', 'full'

    started_at TIMESTAMP,
    completed_at TIMESTAMP,

    papers_processed INTEGER,
    papers_added INTEGER,
    papers_updated INTEGER,
    papers_removed INTEGER,

    errors JSONB,

    created_at TIMESTAMP DEFAULT NOW()
);
```

---

## 7. Resource Estimates

### 7.1 Processing Time Estimates

| Stage | 1M Papers | 5M Papers | 10M Papers |
|-------|-----------|-----------|------------|
| **Download (baseline)** | ~2 hours | ~5 hours | ~10 hours |
| **XML Parsing** | ~20 min | ~1.5 hours | ~3 hours |
| **Filtering** | ~10 min | ~45 min | ~1.5 hours |
| **Entity Extraction** | ~15 min | ~1 hour | ~2 hours |
| **Embedding (CPU)** | ~17 min | ~1.5 hours | ~3 hours |
| **Embedding (GPU)** | ~1 min | ~5 min | ~11 min |
| **RuVector Import** | ~10 min | ~45 min | ~1.5 hours |
| **Supabase Import** | ~15 min | ~1 hour | ~2 hours |
| **Total (CPU)** | **~3 hours** | **~12 hours** | **~24 hours** |
| **Total (GPU)** | **~2.5 hours** | **~10 hours** | **~20 hours** |

### 7.2 Storage Requirements

| Component | 1M Papers | 5M Papers | 10M Papers |
|-----------|-----------|-----------|------------|
| **Raw XML (temp)** | ~3 GB | ~15 GB | ~30 GB |
| **Processed JSON** | ~500 MB | ~2.5 GB | ~5 GB |
| **RuVector (tiered)** | ~400 MB | ~2 GB | ~4 GB |
| **Supabase metadata** | ~500 MB | ~2.5 GB | ~5 GB |
| **HNSW Index** | ~600 MB | ~3 GB | ~6 GB |
| **Total permanent** | **~2 GB** | **~10 GB** | **~20 GB** |

### 7.3 Compute Costs

#### Initial Load (One-time)

| Resource | Duration | Cost |
|----------|----------|------|
| **CPU (t3.medium)** | 24 hours | ~$1.20 |
| **GPU (g4dn.xlarge)** | 2 hours | ~$1.00 |
| **Download bandwidth** | 100 GB | ~$0 (NCBI FTP) |
| **Total** | - | **~$1-2** |

#### Daily Updates (Ongoing)

| Resource | Duration | Cost/Day | Cost/Month |
|----------|----------|----------|------------|
| **CPU processing** | ~5 min | ~$0.01 | ~$0.30 |
| **Supabase writes** | - | ~$0 | ~$0 |
| **Total** | - | **~$0.01** | **~$0.30** |

### 7.4 Complete Cost Summary

| Phase | Database Size | Processing | Storage | Monthly Total |
|-------|---------------|------------|---------|---------------|
| **MVP (1M papers)** | 2 GB | $0.30 | $0 (free tier) | **~$0.30** |
| **Standard (5M)** | 10 GB | $0.30 | $5 | **~$5.30** |
| **Comprehensive (10M)** | 20 GB | $0.30 | $10 | **~$10.30** |

---

## 8. Complete Pipeline Script

```python
#!/usr/bin/env python3
"""
PubMed ETL Pipeline for Gene Knowledge Base

Usage:
    python pipeline.py --mode full      # Full initial load
    python pipeline.py --mode daily     # Daily incremental update
    python pipeline.py --mode weekly    # Weekly maintenance
"""

import argparse
import asyncio
from datetime import datetime
from pathlib import Path
import logging
import os

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Paths
DATA_DIR = Path("/home/ubuntu/src/gene/data")
RAW_DIR = DATA_DIR / "raw/ncbi/pubmed"
PROCESSED_DIR = DATA_DIR / "processed/research"
EMBEDDINGS_DIR = DATA_DIR / "embeddings"

# Ensure directories exist
for dir_path in [RAW_DIR, PROCESSED_DIR, EMBEDDINGS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)


async def run_full_pipeline(max_papers: int = None):
    """
    Run full ETL pipeline.

    1. Download baseline files
    2. Parse and filter
    3. Extract entities
    4. Generate embeddings
    5. Load to RuVector and Supabase
    """
    from extract.pubmed_downloader import PubMedDownloader
    from process.pubmed_parser import PubMedXMLParser
    from filter.genetics_filter import GeneticsFilter
    from process.entity_extractor import BioEntityExtractor
    from embeddings.generator import EmbeddingGenerator
    from embeddings.chunker import ArticleChunker

    logger.info("Starting full ETL pipeline...")
    start_time = datetime.now()

    # Initialize components
    downloader = PubMedDownloader(RAW_DIR)
    parser = PubMedXMLParser()
    genetics_filter = GeneticsFilter()
    entity_extractor = BioEntityExtractor()
    embedding_generator = EmbeddingGenerator()
    chunker = ArticleChunker()

    # Stats
    total_processed = 0
    total_relevant = 0

    # Process each baseline file
    for xml_file in downloader.download_baseline():
        logger.info(f"Processing {xml_file.name}...")

        file_relevant = 0
        batch = []

        for article in parser.parse_file(xml_file):
            total_processed += 1

            # Quick filter first
            if not genetics_filter.quick_filter(f"{article.title} {article.abstract}"):
                continue

            # Full filter
            result = genetics_filter.filter_paper(
                title=article.title,
                abstract=article.abstract,
                mesh_terms=article.mesh_terms
            )

            if result.is_relevant:
                # Extract entities
                entities = entity_extractor.extract_all(
                    f"{article.title} {article.abstract}"
                )

                article.extracted_genes = list(entities.genes)
                article.extracted_snps = list(entities.snps)
                article.relevance_score = result.confidence

                batch.append(article)
                file_relevant += 1
                total_relevant += 1

            # Process in batches
            if len(batch) >= 1000:
                await process_batch(batch, embedding_generator, chunker)
                batch = []

            # Check max papers limit
            if max_papers and total_relevant >= max_papers:
                break

        # Process remaining batch
        if batch:
            await process_batch(batch, embedding_generator, chunker)

        logger.info(f"  - Relevant: {file_relevant}")

        if max_papers and total_relevant >= max_papers:
            logger.info(f"Reached max papers limit: {max_papers}")
            break

    # Final stats
    duration = datetime.now() - start_time
    logger.info(f"Pipeline complete!")
    logger.info(f"  Total processed: {total_processed:,}")
    logger.info(f"  Total relevant: {total_relevant:,}")
    logger.info(f"  Duration: {duration}")


async def process_batch(articles, embedding_generator, chunker):
    """Process a batch of articles."""
    from load.ruvector_loader import RuVectorLoader
    from load.supabase_loader import SupabaseLoader

    # Generate embeddings
    texts = []
    for article in articles:
        chunk = chunker.chunk_for_rag(
            pmid=article.pmid,
            title=article.title,
            abstract=article.abstract
        )
        texts.append(chunk.text)

    embeddings = embedding_generator.generate_batch(texts)

    # Prepare records
    ruvector_records = []
    supabase_records = []

    for article, embedding in zip(articles, embeddings):
        ruvector_records.append({
            'id': f"pmid:{article.pmid}",
            'embedding': embedding.tolist(),
            'metadata': {
                'pmid': article.pmid,
                'title': article.title[:200],
            }
        })

        supabase_records.append({
            'pmid': article.pmid,
            'title': article.title,
            'abstract': article.abstract[:2000],
            'authors': article.authors,
            'journal': article.journal,
            'publication_date': article.publication_date,
            'mesh_terms': article.mesh_terms,
            'extracted_genes': article.extracted_genes,
            'extracted_snps': article.extracted_snps,
            'relevance_score': article.relevance_score,
        })

    # Load to databases
    # RuVectorLoader().import_batch(ruvector_records)
    # SupabaseLoader().upsert_articles(supabase_records)


async def run_daily_update():
    """Run daily incremental update."""
    from updates.update_manager import UpdateManager

    manager = UpdateManager(DATA_DIR / "state")
    result = manager.run_daily_update()

    logger.info(f"Daily update complete: {result}")


async def run_weekly_maintenance():
    """Run weekly maintenance tasks."""
    from updates.update_manager import UpdateManager

    manager = UpdateManager(DATA_DIR / "state")
    result = manager.run_weekly_check()

    logger.info(f"Weekly maintenance complete: {result}")


def main():
    parser = argparse.ArgumentParser(description='PubMed ETL Pipeline')
    parser.add_argument('--mode', choices=['full', 'daily', 'weekly'],
                        required=True, help='Pipeline mode')
    parser.add_argument('--max-papers', type=int, default=None,
                        help='Max papers to process (for testing)')

    args = parser.parse_args()

    if args.mode == 'full':
        asyncio.run(run_full_pipeline(args.max_papers))
    elif args.mode == 'daily':
        asyncio.run(run_daily_update())
    elif args.mode == 'weekly':
        asyncio.run(run_weekly_maintenance())


if __name__ == '__main__':
    main()
```

---

## 9. Cron Schedule

```bash
# /etc/cron.d/gene-papers-etl

# Daily updates (2 AM, after PubMed daily release)
0 2 * * * ubuntu /usr/bin/python3 /home/ubuntu/src/gene/scripts/etl/papers/pipeline.py --mode daily >> /var/log/gene-etl/daily.log 2>&1

# Weekly maintenance (Sunday 3 AM)
0 3 * * 0 ubuntu /usr/bin/python3 /home/ubuntu/src/gene/scripts/etl/papers/pipeline.py --mode weekly >> /var/log/gene-etl/weekly.log 2>&1

# Monthly full refresh check (1st of month, 4 AM)
0 4 1 * * ubuntu /usr/bin/python3 /home/ubuntu/src/gene/scripts/etl/papers/check_refresh.py >> /var/log/gene-etl/monthly.log 2>&1
```

---

## 10. Monitoring and Alerts

```python
# scripts/etl/monitoring/health_check.py

from datetime import datetime, timedelta
import os

def check_pipeline_health() -> dict:
    """Check ETL pipeline health status."""
    from supabase import create_client

    client = create_client(os.environ['SUPABASE_URL'], os.environ['SUPABASE_KEY'])

    # Check last update
    last_update = client.table('etl_update_log') \
        .select('completed_at') \
        .order('completed_at', desc=True) \
        .limit(1) \
        .execute()

    last_update_time = datetime.fromisoformat(last_update.data[0]['completed_at'])
    hours_since_update = (datetime.now() - last_update_time).total_seconds() / 3600

    # Check paper count
    paper_count = client.table('research_articles') \
        .select('pmid', count='exact') \
        .execute()

    return {
        'status': 'healthy' if hours_since_update < 48 else 'stale',
        'last_update': last_update_time.isoformat(),
        'hours_since_update': hours_since_update,
        'total_papers': paper_count.count,
        'alerts': [] if hours_since_update < 48 else ['No update in 48 hours']
    }
```

---

## Summary

| Component | Choice | Rationale |
|-----------|--------|-----------|
| **Source** | PubMed FTP baseline + daily | Free, comprehensive, well-structured |
| **Parsing** | lxml (streaming) | 10x faster than BeautifulSoup, memory-efficient |
| **Filtering** | MeSH + keyword + entity | 85% reduction to relevant papers |
| **Embedding model** | all-MiniLM-L6-v2 | RuVector optimized, 384d, fast |
| **Vector storage** | RuVector | 75% compression, graph queries |
| **Metadata storage** | Supabase PostgreSQL | Managed, generous free tier |
| **Update frequency** | Daily | Balance freshness vs cost |
| **Initial load** | ~24 hours CPU / ~20 hours GPU | One-time cost |
| **Total storage** | 2-20 GB (RuVector compressed) | Depends on paper count |
| **Monthly cost** | ~$0.30-10 | Mostly free with Supabase free tier |
