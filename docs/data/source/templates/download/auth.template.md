## Authentication

### Registration

1. Go to {{REGISTRATION_URL}}
2. Create an account with {{REQUIREMENTS}}
3. Verify your email
4. Navigate to {{API_KEY_PATH}}
5. Generate API key

### API Key Management

| Property | Value |
|----------|-------|
| Key Format | `{{KEY_FORMAT}}` |
| Expiration | {{EXPIRATION}} |
| Regeneration | {{Yes/No}} |

### Usage

#### Header Authentication

```bash
curl -H "Authorization: Bearer {{API_KEY}}" {{URL}}
```

#### Query Parameter

```bash
curl "{{URL}}?api_key={{API_KEY}}"
```

### Environment Variable

```bash
export {{SOURCE_NAME}}_API_KEY="your-key-here"
```

### Rate Limits by Auth Level

| Level | Requests/hour | Requests/day |
|-------|---------------|--------------|
| Anonymous | {{N}} | {{N}} |
| Basic | {{N}} | {{N}} |
| Premium | {{N}} | {{N}} |

### Troubleshooting

| Error | Cause | Solution |
|-------|-------|----------|
| 401 Unauthorized | Invalid key | Regenerate API key |
| 403 Forbidden | Insufficient permissions | Upgrade account |
| 429 Too Many Requests | Rate limited | Wait or upgrade |
