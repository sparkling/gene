---
id: download-natural-medicines
title: "Natural Medicines Download Instructions"
type: download
parent: _index.md
last_updated: 2026-01-24
access_type: restricted
---

# Natural Medicines Download Instructions

## Access Classification

| Aspect | Status |
|--------|--------|
| Access Type | **Restricted - Subscription Required** |
| Bulk Downloads | Not Available |
| API Access | Enterprise License Only |
| Data Export | Limited (Patient handouts, PDFs) |
| EHR Integration | Available (Epic, Cerner, etc.) |
| Redistribution | Prohibited |

---

## Important Notice

Natural Medicines is a proprietary subscription database. **Bulk data downloads are not available.** Data must be accessed through the web interface or EHR integration with a valid institutional or individual subscription.

## Access Requirements

| Requirement | Details |
|-------------|---------|
| Subscription | Required (institutional or individual) |
| Account | Registration required |
| Terms of Use | No bulk downloading, no redistribution |

## Subscription Options

| Plan | Pricing | Features |
|------|---------|----------|
| Individual | Contact for pricing | Full web access, mobile app |
| Institutional | Contact for pricing | Multi-user, IP authentication |
| Enterprise | Contact for pricing | API access, EHR integration |

Contact: https://trchealthcare.com/

## Available Data Access Methods

### Method 1: Web Interface (Primary)

```
1. Navigate to: https://naturalmedicines.therapeuticresearch.com/
2. Log in with subscription credentials
3. Search for ingredients, interactions, or conditions
4. View monographs and use interaction checker
```

### Method 2: Mobile App

```
1. Download Natural Medicines app (iOS/Android)
2. Log in with subscription credentials
3. Search and browse on mobile device
4. Access offline monographs (limited)
```

### Method 3: EHR Integration (Enterprise)

```
Supported systems:
- Epic
- Cerner
- Allscripts
- Other EHR vendors

Contact TRC Healthcare for integration details.
```

### Method 4: Patient Handouts Export

```
1. Open monograph in web interface
2. Click "Patient Handout" option
3. Select language and detail level
4. Download PDF or print
```

## Data Export Options

### Available Exports (Limited)

| Data Type | Export Method | Format |
|-----------|---------------|--------|
| Patient handouts | Download | PDF |
| Interaction reports | Print/Save | PDF |
| Individual monographs | Print | PDF |

### Not Available

- Bulk database downloads
- API access (without enterprise license)
- Raw data exports
- Automated data extraction

## Alternative Data Sources

For programmatic access to natural product data, consider:

| Source | Access | Coverage |
|--------|--------|----------|
| [DSLD](../dsld/) | Free API | US supplement labels |
| [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | Free | Compound data |
| [DrugBank](https://go.drugbank.com/) | Free/Paid | Drug interactions |
| [HMDB](../../6.4.metabolomics/hmdb/) | Free | Metabolites |

## Terms of Service Highlights

- Subscription required for access
- No systematic downloading
- No redistribution
- Citation permitted with attribution
- Patient handouts may be shared with patients

## Manual Research Workflow

```bash
# For personal research notes (comply with ToS)

# 1. Create structured notes
mkdir -p nm_research/{ingredients,interactions,conditions}

# 2. Document findings manually
cat > nm_research/ingredients/turmeric_notes.md << 'EOF'
# Turmeric Research Notes
Source: Natural Medicines Database
Access Date: 2026-01-23
Subscription Required: Yes

## Key Information
- Safety Grade: Likely Safe (oral, appropriate doses)
- Effectiveness: Possibly Effective for osteoarthritis

## Notable Drug Interactions
- Anticoagulants: May increase bleeding risk
- Antidiabetic drugs: May enhance hypoglycemic effect

## Citation
Natural Medicines. Turmeric.
Natural Medicines Comprehensive Database.
TRC Healthcare. [Access date]
EOF
```

## Interaction Checker Usage

```
1. Navigate to Interaction Checker tool
2. Enter patient's medications
3. Enter supplements/natural products
4. Review severity ratings and recommendations
5. Print/save report for clinical use
```

## Evidence Rating Reference

| Rating | Use in Research |
|--------|-----------------|
| Effective | Strong recommendation basis |
| Likely Effective | Good recommendation basis |
| Possibly Effective | May recommend with caveats |
| Insufficient Evidence | Research opportunity |

## Update Schedule

| Content Type | Update Frequency |
|--------------|------------------|
| Monograph updates | Continuous (as evidence emerges) |
| New monographs | As needed |
| Interaction database | Continuous |
| Safety alerts | Immediate |

## Institutional Access Options

### IP Authentication
```
1. Contact TRC Healthcare
2. Provide institutional IP range
3. Users access without individual login from campus
```

### Proxy Authentication
```
1. Configure proxy server
2. Integrate with library authentication
3. Off-campus users authenticate via proxy
```

### SAML/SSO Integration
```
1. Configure identity provider
2. Integrate with institutional SSO
3. Seamless authentication for users
```

## Common Questions

### Can I download the database?
No, Natural Medicines does not provide bulk downloads. Data is proprietary.

### Is there a free version?
No, subscription is required for access. Some institutions provide access to students/faculty.

### Can I use data in publications?
Yes, with proper citation. Contact publisher for extensive use.

### What about API access?
API access requires enterprise license. Contact TRC Healthcare sales.

## Support Contacts

| Type | Contact |
|------|---------|
| Technical Support | support@trchealthcare.com |
| Subscription Sales | sales@trchealthcare.com |
| Institutional | institutions@trchealthcare.com |

---

## Data Structure Overview (For Subscribers)

### Available Data Categories

| Category | Content | Export Options |
|----------|---------|----------------|
| Ingredient Monographs | 1,400+ natural products | View/Print/PDF |
| Drug Interactions | 175,000+ interaction pairs | Interaction Checker reports |
| Effectiveness Ratings | Evidence-based ratings by condition | View/Print |
| Patient Handouts | Consumer-friendly summaries | PDF (customizable) |
| Continuing Education | CE/CME modules | Online only |

### Data Fields Available Through Interface

**Monograph Fields:**
- Ingredient name (common, scientific, synonyms)
- Taxonomy and botanical family
- Parts used, traditional uses
- Effectiveness ratings by condition (Effective to Insufficient Evidence)
- Safety ratings (Likely Safe to Unsafe)
- Mechanism of action
- Pharmacokinetics (absorption, metabolism, half-life)
- Adverse effects and contraindications
- Dosing recommendations by condition
- Primary literature references

**Interaction Fields:**
- Natural product and drug names
- Severity rating (Contraindicated, Major, Moderate, Minor, Unknown)
- Documentation level (Well-documented, Theoretical, etc.)
- Mechanism of interaction
- Clinical significance and management recommendations
- Supporting references (PubMed IDs)

---

## Programmatic Alternatives

For researchers needing programmatic access to natural product evidence data, consider these open alternatives:

### DrugBank (Drug Interactions)
```bash
# Free tier available for academic use
# Natural product interactions available
curl "https://go.drugbank.com/unearth/q?searcher=drugs&query=curcumin" \
  -H "Accept: application/json"
```

### PubChem BioAssay
```bash
# Free access to bioactivity data
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/curcumin/assaysummary/JSON"
```

### Open Targets Platform
```bash
# Free drug-target-disease associations
curl "https://api.platform.opentargets.org/api/v4/graphql" \
  -X POST \
  -H "Content-Type: application/json" \
  -d '{"query": "query { search(queryString: \"curcumin\") { total hits { id name } } }"}'
```

### NIH ODS Fact Sheets
```bash
# Free supplement information
# Web: https://ods.od.nih.gov/factsheets/list-all/
# No API, but structured HTML can be parsed
```

---

## Related Resources

- [Schema Documentation](./schema.md)
- [DSLD](../dsld/) - Free government supplement database
- [ConsumerLab](../consumerlab/) - Product testing database
