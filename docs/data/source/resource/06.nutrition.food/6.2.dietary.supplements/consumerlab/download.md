---
id: download-consumerlab
title: "ConsumerLab Download Instructions"
type: download
parent: README.md
last_updated: 2026-01-24
access_type: restricted
---

# ConsumerLab Download Instructions

## Access Classification

| Aspect | Status |
|--------|--------|
| Access Type | **Restricted - Subscription Required** |
| Bulk Downloads | Not Available |
| API Access | Not Public |
| Data Export | Limited (PDF/Print only) |
| Redistribution | Prohibited |

---

## Important Notice

ConsumerLab is a proprietary subscription service. **Bulk data downloads are not available.** Data must be accessed through the web interface with a valid subscription.

## Access Requirements

| Requirement | Details |
|-------------|---------|
| Subscription | $49.95/year (individual) |
| Account | Email registration required |
| Terms of Use | No bulk downloading, no redistribution |

## Subscription Options

| Plan | Price | Features |
|------|-------|----------|
| Individual | $49.95/year | Full access to reviews and monographs |
| Institutional | Contact sales | Multi-user access, API inquiries |

## Available Data Access Methods

### Method 1: Web Interface (Primary)

```
1. Navigate to: https://www.consumerlab.com
2. Log in with subscription credentials
3. Browse or search for products
4. View individual test results and reviews
```

### Method 2: Manual Data Collection

```bash
# Save individual pages (for personal use only)
# Requires active login session

# Example: Save a product review page
curl -b cookies.txt "https://www.consumerlab.com/reviews/vitamin-d/" -o vitamin_d_review.html

# Note: Check Terms of Service before automated access
```

### Method 3: Contact for Institutional Access

```
For bulk data licensing:
- Email: info@consumerlab.com
- Phone: 914-722-9149
- Web form: https://www.consumerlab.com/contact/
```

## Data Export Options

### Available Exports (Manual)

| Data Type | Export Method | Format |
|-----------|---------------|--------|
| Product comparisons | Print/PDF | PDF |
| Review summaries | Copy/paste | Text |
| Encyclopedia entries | Print | PDF |

### Not Available

- Bulk data downloads
- API access (public)
- Database exports
- Automated scraping (prohibited)

## Alternative Data Sources

For programmatic access to supplement data, consider:

| Source | Access | Coverage |
|--------|--------|----------|
| [DSLD](../dsld/) | Free API | US supplement labels |
| [NIH ODS](https://ods.od.nih.gov/) | Free | Supplement fact sheets |
| [Examine.com](https://examine.com/) | Free/Paid | Evidence summaries |

## Terms of Service Highlights

- Personal use only
- No redistribution
- No systematic downloading
- No data mining
- Citation permitted with attribution

## Manual Data Collection Workflow

```bash
# For personal research notes (comply with ToS)

# 1. Create local notes structure
mkdir -p consumerlab_notes/{vitamins,minerals,herbs}

# 2. Record key findings manually
cat > consumerlab_notes/vitamins/vitamin_d_summary.md << 'EOF'
# Vitamin D Product Summary
Source: ConsumerLab.com
Access Date: 2026-01-23
Subscription Required: Yes

## Key Findings
- [Note your findings here]

## Top Approved Products
- [List products manually]

## Citation
ConsumerLab.com. Vitamin D Supplements Review.
https://www.consumerlab.com/reviews/vitamin-d/
Accessed: 2026-01-23
EOF
```

## Update Schedule

| Content Type | Update Frequency |
|--------------|------------------|
| Product reviews | As tested (ongoing) |
| Encyclopedia | Continuously updated |
| Drug interactions | As evidence emerges |
| Quality alerts | Immediate |

## Common Questions

### Can I download the database?
No, ConsumerLab does not provide bulk data downloads. Data is proprietary and accessible only through subscription.

### Is there an API?
No public API is available. Institutional customers may inquire about custom data access.

### Can I use data in publications?
Yes, with proper citation. Contact ConsumerLab for commercial use licensing.

### What alternatives exist for bulk data?
- DSLD (NIH) - Free API for US supplement labels
- PubMed - Literature on supplement research
- NIH ODS - Factsheets on dietary supplements

---

## Data Structure Overview (For Subscribers)

### Available Data Categories

| Category | Content | Export Options |
|----------|---------|----------------|
| Product Reviews | Test results, approval status, cost analysis | View/Print/PDF |
| Ingredient Encyclopedia | Evidence summaries, dosing, interactions | View/Print/PDF |
| Brand Comparisons | Side-by-side test results | View/Print |
| Quality Alerts | Recalls, warnings, contamination notices | View only |
| Drug Interactions | Supplement-drug interaction checker | View/Print |

### Data Fields Available Through Interface

**Product Review Fields:**
- Product name, brand, UPC
- Approval status (Approved/Not Approved)
- Test date and results
- Potency (% of label claim)
- Heavy metal content (lead, arsenic, cadmium, mercury)
- Disintegration time
- Contaminant screening results
- Cost per serving analysis

**Ingredient Monograph Fields:**
- Common/scientific names
- Evidence ratings by health claim
- Recommended intake ranges
- Drug interaction warnings
- Form comparisons (e.g., D2 vs D3)
- Safety considerations

---

## Programmatic Alternatives

For researchers needing programmatic access to supplement quality data, consider these open alternatives:

### DSLD (Dietary Supplement Label Database)
```bash
# Free API access to US supplement label data
curl "https://api.ods.od.nih.gov/dsld/v9/browse/ingredient?name=vitamin%20d" \
  -H "Accept: application/json"

# Includes: product names, ingredients, amounts, manufacturers
```

### NSF Certified Products Database
```bash
# Search NSF-certified supplements (free)
# Web interface: https://www.nsf.org/certified-products-systems
# Includes: certified products, testing standards
```

### USP Verified Products
```bash
# Search USP-verified supplements (free)
# Web interface: https://www.quality-supplements.org/verified-products
# Includes: verified products, quality standards
```

---

## Related Resources

- [Schema Documentation](./schema.md)
- [DSLD](../dsld/) - Free government supplement database
- [Natural Medicines](../natural.medicines/) - Evidence-based database
