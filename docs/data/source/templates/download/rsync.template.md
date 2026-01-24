## Rsync Mirroring

### Rsync URL

```
rsync://{{RSYNC_HOST}}/{{MODULE}}/
```

### Download Commands

#### List Contents

```bash
rsync --list-only rsync://{{RSYNC_HOST}}/{{MODULE}}/
```

#### Sync Directory

```bash
rsync -avz --progress rsync://{{RSYNC_HOST}}/{{MODULE}}/ ./local/
```

#### Incremental Update

```bash
rsync -avz --delete rsync://{{RSYNC_HOST}}/{{MODULE}}/ ./local/
```

### Bandwidth Limiting

```bash
rsync -avz --bwlimit={{KBPS}} rsync://{{RSYNC_HOST}}/{{MODULE}}/ ./local/
```

### Mirror Sites

| Location | URL |
|----------|-----|
| {{LOCATION_1}} | `rsync://{{HOST_1}}/{{MODULE}}/` |
| {{LOCATION_2}} | `rsync://{{HOST_2}}/{{MODULE}}/` |

### Recommended Cron Schedule

```bash
# Sync weekly on Sunday at 2am
0 2 * * 0 rsync -avz --delete rsync://{{RSYNC_HOST}}/{{MODULE}}/ /data/{{SOURCE}}/
```
