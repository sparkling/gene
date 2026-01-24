# {{SOURCE_NAME}} Download & Access

## Quick Start

```bash
{{QUICK_DOWNLOAD_COMMAND}}
```

## Access Methods

| Method | Available | Details |
|--------|-----------|---------|
| Direct Download | {{Yes/No}} | {{URL or N/A}} |
| REST API | {{Yes/No}} | See [REST API](#rest-api) |
| SPARQL | {{Yes/No}} | See [SPARQL](#sparql) |
| FTP/Bulk | {{Yes/No}} | See [Bulk Download](#bulk-download) |

## Available Formats

| Format | Extension | Compression | Description |
|--------|-----------|-------------|-------------|
| {{FORMAT_1}} | `.{{ext}}` | {{gzip/none}} | {{Description}} |
| {{FORMAT_2}} | `.{{ext}}` | {{gzip/none}} | {{Description}} |

---

## Dataset Versions

### Current Release: {{VERSION}}

| Property | Value |
|----------|-------|
| Version | {{VERSION}} |
| Release Date | {{YYYY-MM-DD}} |
| Total Size | {{SIZE}} |
| Checksum | `sha256:{{HASH}}` |

### Version Contents

| Component | Size | Records | Description |
|-----------|------|---------|-------------|
| {{FILE_1}} | {{SIZE}} | {{COUNT}} | {{Description}} |
| {{FILE_2}} | {{SIZE}} | {{COUNT}} | {{Description}} |

### Previous Versions

| Version | Release | Size | Status |
|---------|---------|------|--------|
| {{PREV_VERSION}} | {{DATE}} | {{SIZE}} | Archived |

---

## Update Schedule

| Property | Value |
|----------|-------|
| Frequency | {{FREQUENCY}} |
| Release Day | {{DAY}} |
| Notification | {{URL/RSS}} |

---

<!-- Include applicable sections below based on Access Methods table -->

<!--
## REST API
(Include from download/rest-api.template.md if applicable)
-->

<!--
## SPARQL
(Include from download/sparql.template.md if applicable)
-->

<!--
## Bulk Download
(Include from download/ftp.template.md if applicable)
-->

---

## Post-Download Processing

```bash
# Decompress (if needed)
{{DECOMPRESS_COMMAND}}

# Validate checksum
{{CHECKSUM_COMMAND}}

# Index (if applicable)
{{INDEX_COMMAND}}
```

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| {{COMMON_ISSUE_1}} | {{SOLUTION_1}} |
| {{COMMON_ISSUE_2}} | {{SOLUTION_2}} |

---

## See Also

- [README](./README.md) - Source overview
- [Schema](./schema.json) - Data structure
- [Sample](./sample.json) - Example records
