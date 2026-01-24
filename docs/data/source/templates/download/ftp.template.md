## Bulk Download (FTP)

### FTP URL

```
ftp://{{FTP_HOST}}/{{PATH}}/
```

### Directory Structure

```
{{FTP_HOST}}/
├── {{DIR_1}}/
│   ├── {{FILE_1}}
│   └── {{FILE_2}}
└── {{DIR_2}}/
    └── {{FILE_3}}
```

### Download Commands

#### Single File

```bash
wget ftp://{{FTP_HOST}}/{{PATH}}/{{FILE}}
```

#### Recursive Download

```bash
wget -r -np -nH --cut-dirs={{N}} ftp://{{FTP_HOST}}/{{PATH}}/
```

#### With lftp (resumable)

```bash
lftp -c "mirror --parallel=4 ftp://{{FTP_HOST}}/{{PATH}}/ ./local/"
```

### Mirror Sites

| Location | URL | Sync Frequency |
|----------|-----|----------------|
| {{LOCATION_1}} | {{URL_1}} | {{FREQ}} |
| {{LOCATION_2}} | {{URL_2}} | {{FREQ}} |

### Checksums

```bash
# Download and verify
wget ftp://{{FTP_HOST}}/{{PATH}}/{{CHECKSUM_FILE}}
sha256sum -c {{CHECKSUM_FILE}}
```
