# Web Data Sources & Scraping Targets for Gene Platform

**Last Updated:** January 2026
**Purpose:** Comprehensive analysis of web scraping targets and data sources for health/genetics community data, supplement reviews, clinical trials, patient-reported outcomes, health apps, and wearable data.

---

## Executive Summary

| Category | Top Sources | Feasibility | Legal Risk |
|----------|-------------|-------------|------------|
| Health Forums | Reddit (r/genetics, r/23andme), HealthUnlocked | Medium-High | Medium-High |
| Supplement Reviews | ConsumerLab, Labdoor, Examine.com | Low | High |
| Clinical Trials | ClinicalTrials.gov, AACT, WHO ICTRP | High | Low |
| Patient Outcomes | PROMIS, PatientsLikeMe, Open Humans | Medium | Medium |
| Health Apps | Apple HealthKit, Health Connect | Medium | Medium |
| Wearable Data | Open Wearables, Terra API, Fitabase | High | Low-Medium |

---

## 1. Health Forums with Genetic Discussions

### 1.1 Reddit Genetics Communities

| Subreddit | Members | Content Type | Focus |
|-----------|---------|--------------|-------|
| r/23andme | 81,400+ | User results, interpretations | DTC genetic testing |
| r/AncestryDNA | ~50,000 | Ancestry results, ethnicity | Genealogy genetics |
| r/genetics | ~100,000 | Research, testing, pharmacogenomics | Scientific discussion |
| r/genomics | ~30,000 | Research breakthroughs, tools | Technical genomics |
| r/GeneticTesting | ~15,000 | Testing methods, results | Testing guidance |
| r/DNA | ~20,000 | General DNA topics | Broad genetics |

**URL:** https://reddit.com/r/{subreddit}

**Content Type:**
- User-submitted genetic test results and interpretations
- Discussions about SNP meanings and health implications
- Pharmacogenomics questions and experiences
- Treatment decisions based on genetic data

**Data Structure:**
```json
{
  "post": {
    "title": "string",
    "body": "string (markdown)",
    "author": "string",
    "created_utc": "timestamp",
    "score": "integer",
    "comments": [{
      "body": "string",
      "author": "string",
      "score": "integer"
    }]
  }
}
```

**Scraping Feasibility:** MEDIUM
- Official API available via [r/reddit4researchers](https://reddit.com/r/reddit4researchers)
- Research API access requires proposal approval
- Rate limits: 60 requests/minute
- Historical data via Pushshift (limited availability)

**Legal Considerations:** HIGH RISK
- Reddit Terms of Service prohibit unauthorized scraping
- Reddit sued Perplexity AI (October 2025) for DMCA circumvention
- Legitimate research access requires formal application
- Non-commercial academic use may be permitted with approval
- Must anonymize data in publications and credit Reddit

**Recommendation:**
- Apply via r/reddit4researchers for research API access
- Use Reddit's official research partnership program
- Avoid third-party scrapers (SerpApi, Oxylabs named in lawsuits)

**References:**
- [Reddit Developer Platform](https://support.reddithelp.com/hc/en-us/articles/14945211791892-Developer-Platform-Accessing-Reddit-Data)
- [Responsible Builder Policy](https://support.reddithelp.com/hc/en-us/articles/42728983564564-Responsible-Builder-Policy)
- [Genetic Life Hacks - Reddit Forums on Genetics](https://www.geneticlifehacks.com/reddit-forums-on-genetics/)

---

### 1.2 HealthUnlocked

**URL:** https://healthunlocked.com/

**Content Type:**
- Patient support communities (~300+ health conditions)
- Symptom discussions and treatment experiences
- Medication side effects and efficacy reports

**Data Structure:**
- Community-based organization by condition
- Posts with replies, likes, and tags
- User profiles with condition information

**Scraping Feasibility:** LOW
- No public API available
- Platform designed for patient privacy
- Content behind authentication

**Legal Considerations:** HIGH RISK
- Medical data privacy concerns (HIPAA-adjacent)
- Terms of service prohibit data collection
- Patient expectation of privacy

**Recommendation:** NOT RECOMMENDED for scraping. Consider partnership approach.

**Reference:** [HealthUnlocked Communities](https://healthunlocked.com/communities)

---

### 1.3 MyFitnessPal Community

**URL:** https://community.myfitnesspal.com/

**Content Type:**
- Nutrition and diet discussions
- Fitness tracking experiences
- Weight management support

**Data Structure:**
- Forum-based with categories and groups
- User discussions with replies
- Recent integration with ChatGPT Health (January 2026)

**Scraping Feasibility:** LOW
- No public API for community data
- Login required for most content

**Legal Considerations:** MEDIUM RISK
- Under Armour/MyFitnessPal privacy policies
- User data expectations of privacy

**Reference:** [MyFitnessPal Community](https://community.myfitnesspal.com/en)

---

## 2. Supplement Review Sites with Structured Data

### 2.1 ConsumerLab

**URL:** https://www.consumerlab.com/

**Content Type:**
- Independent lab testing results
- Product purity analysis
- Label accuracy verification
- Contaminant testing

**Data Structure:**
```
Product Reviews:
- Product name, brand, category
- Label accuracy score
- Purity test results
- Contaminant levels (when tested)
- Overall rating
- Price/value assessment
```

**Scraping Feasibility:** LOW
- Paid subscription required ($99.95/2 years)
- No public API
- Content behind paywall
- Anti-scraping measures likely

**Legal Considerations:** HIGH RISK
- Copyrighted test results
- Commercial product (subscription model)
- Terms of service restrictions

**Recommendation:** NOT RECOMMENDED for scraping. Purchase subscription for internal research use only.

**Reference:** [ConsumerLab Reviews](https://www.consumerlab.com/reviews/)

---

### 2.2 Labdoor

**URL:** https://labdoor.com/

**Content Type:**
- Third-party lab testing
- Quality scores (Label Accuracy, Purity, Nutritional Value, Safety, Efficacy)
- Price rankings and value scores

**Data Structure:**
```
Product Ratings:
- Overall Quality Score (0-100)
- Label Accuracy Score
- Product Purity Score
- Nutritional Value Score
- Ingredient Safety Score
- Projected Efficacy Score
- Price comparison
```

**Scraping Feasibility:** MEDIUM
- Public rankings visible without login
- No official API available
- Structured HTML data

**Legal Considerations:** HIGH RISK
- Proprietary scoring methodology
- Terms of service restrictions
- Commercial entity (VC-backed: Rock Health, Mark Cuban, Y Combinator)

**Recommendation:** Monitor for API availability. Do not scrape without permission.

**Reference:** [Labdoor About](https://labdoor.com/about)

---

### 2.3 Examine.com

**URL:** https://examine.com/

**Content Type:**
- Evidence-based supplement analysis
- Research summaries by compound
- Dosage recommendations
- Effect ratings with confidence levels

**Data Structure:**
```
Supplement Pages:
- Compound name and aliases
- Health effects (rated by evidence level)
- Human studies summary
- Dosage recommendations
- Interaction warnings
- Research citations
```

**Scraping Feasibility:** LOW-MEDIUM
- Subscription required for full access
- Public summaries available
- No official API

**Legal Considerations:** HIGH RISK
- Copyrighted editorial content
- Subscription-based business model
- Attribution requirements

**Recommendation:** Use PubMed/NIH ODS as alternative for similar evidence-based data. See existing `data-sources-public.md` for NIH ODS API details.

**Reference:** [Examine.com](https://examine.com/)

---

## 3. Clinical Trial Registries

### 3.1 ClinicalTrials.gov (PRIMARY)

**URL:** https://clinicaltrials.gov/

**API URL:** https://clinicaltrials.gov/api/

**Content Type:**
- Protocol information (study design, eligibility, interventions)
- Results data (outcomes, adverse events)
- Study status and recruitment info
- Sponsor and investigator details

**Data Structure (API v2.0):**
```json
{
  "protocolSection": {
    "identificationModule": {
      "nctId": "NCT00000000",
      "briefTitle": "string",
      "officialTitle": "string"
    },
    "statusModule": {
      "overallStatus": "RECRUITING|COMPLETED|...",
      "startDate": "ISO 8601",
      "completionDate": "ISO 8601"
    },
    "conditionsModule": {
      "conditions": ["string"]
    },
    "interventionsModule": {
      "interventions": [{
        "type": "DRUG|DEVICE|...",
        "name": "string",
        "description": "string"
      }]
    },
    "outcomesModule": {
      "primaryOutcomes": [{}],
      "secondaryOutcomes": [{}]
    },
    "eligibilityModule": {
      "criteria": "string",
      "healthyVolunteers": "boolean",
      "sex": "ALL|FEMALE|MALE",
      "minimumAge": "string",
      "maximumAge": "string"
    }
  },
  "resultsSection": {
    "participantFlowModule": {},
    "baselineCharacteristicsModule": {},
    "outcomeMeasuresModule": {},
    "adverseEventsModule": {}
  }
}
```

**API Features (v2.0 - Current):**
- REST API with OpenAPI 3.0 specification
- JSON primary response format
- ISO 8601 dates
- CommonMark markdown for rich text
- Enumerated values for key fields
- Classic API retired June 2024

**Scraping Feasibility:** HIGH
- Official REST API freely available
- No registration required for basic access
- Bulk download available
- Daily updates

**Legal Considerations:** LOW RISK
- Public domain data
- NIH/NLM operated
- Open access encouraged
- Must comply with Terms and Conditions

**Tools:**
- [ctrdata R Package](https://cran.r-project.org/package=ctrdata) - v1.26.0.9000 (reviewed Jan 2026)
- [clinicaltrialr R Package](https://github.com/serghiou/clinicaltrialr)
- MCP Server Tool available

**References:**
- [ClinicalTrials.gov API Documentation](https://clinicaltrials.gov/data-api/api)
- [API v2.0 Announcement](https://www.nlm.nih.gov/pubs/techbull/ma24/ma24_clinicaltrials_api.html)

---

### 3.2 AACT Database (Recommended for Bulk Access)

**URL:** https://aact.ctti-clinicaltrials.org/

**Content Type:**
- Complete relational database of ClinicalTrials.gov
- Daily updates from source
- PostgreSQL format

**Data Structure:**
- 50+ tables covering all study data elements
- Normalized relational structure
- Pre-processed for analysis

**Download Options:**
- Direct cloud database access (PostgreSQL)
- Daily pg_dump snapshots
- Connection limit: 10 concurrent per account

**Scraping Feasibility:** HIGH (Direct Download)
- Free account registration
- Cloud-based access
- Downloadable snapshots: https://aact.ctti-clinicaltrials.org/snapshots

**Legal Considerations:** LOW RISK
- Public domain derivative
- Free for research use
- Maintained by Clinical Trials Transformation Initiative (CTTI)

**Recommendation:** Use AACT for bulk analysis. Provides cleaner relational structure than raw API.

**Reference:** [AACT Database GitHub](https://github.com/ctti-clinicaltrials/aact)

---

### 3.3 EU Clinical Trials Register (CTIS/EudraCT)

**URLs:**
- CTIS: https://euclinicaltrials.eu/
- EudraCT: https://eudract.ema.europa.eu/
- Search: https://www.clinicaltrialsregister.eu/ctr-search/search

**Content Type:**
- EU/EEA clinical trial data
- Protocol information
- Results data (post January 2022)

**Data Structure:**
- Similar to ClinicalTrials.gov
- EudraCT number identifiers
- Multi-language support

**Current Status (2026):**
- CTIS became single entry point January 31, 2025
- All ongoing trials transitioned from EudraCT
- CTIS designated as WHO ICTRP primary registry
- Clinical trial map launched March 2025
- Task force implementing proposed changes from 2026

**Scraping Feasibility:** MEDIUM
- Public search portal available
- No documented public API
- Web scraping technically possible but challenging

**Legal Considerations:** LOW-MEDIUM RISK
- Public data under EU transparency regulations
- Clinical Trials Regulation requires public availability
- Some data may be exempted for commercial interests

**Recommendation:** Use search portal for targeted queries. Monitor for API development.

**References:**
- [CTIS EMA Page](https://www.ema.europa.eu/en/human-regulatory-overview/research-development/clinical-trials-human-medicines/clinical-trials-information-system)
- [EU Clinical Trials Register](https://www.clinicaltrialsregister.eu/)

---

### 3.4 WHO ICTRP (International Clinical Trials Registry Platform)

**URL:** https://www.who.int/tools/clinical-trials-registry-platform

**Search Portal:** https://trialsearch.who.int

**Content Type:**
- Aggregated data from 17+ national registries
- Global clinical trial coverage
- WHO 20-item Trial Registration Data Set

**Data Structure:**
- Standardized across source registries
- Trial ID, Title, Intervention, Conditions
- Recruitment status and dates
- Primary and secondary outcomes

**Data Access:**
- Weekly updates
- COVID-19 subset: Monthly download files
- Full dataset: Request via Sharepoint form
- Crawling service: Currently unavailable (survey for 2025+ access)

**Scraping Feasibility:** MEDIUM
- Public search portal
- Bulk download via request form
- No REST API currently

**Legal Considerations:** LOW RISK
- WHO public health mandate
- Attribution required
- Non-commercial use clause

**Terms of Use:**
1. Attribute source as WHO ICTRP
2. Keep data current
3. Display processing date
4. No marketing/commercial use

**Reference:** [ICTRP Data Download](https://www.who.int/tools/clinical-trials-registry-platform/network/who-data-set/downloading-records-from-the-ictrp-database)

---

## 4. Patient-Reported Outcome Databases

### 4.1 PROMIS (Patient-Reported Outcomes Measurement Information System)

**URL:** https://www.healthmeasures.net/explore-measurement-systems/promis

**Content Type:**
- Standardized PRO measures
- Physical, mental, social health domains
- Computer adaptive testing (CAT) instruments
- Adult and pediatric versions

**Domains Covered:**
- Pain
- Fatigue
- Emotional distress
- Physical functioning
- Social role participation
- Sleep disturbance
- Anxiety/Depression

**Data Structure:**
- Item banks with calibrated questions
- Standardized T-scores (mean=50, SD=10)
- Domain-specific short forms
- Full item bank parameters available

**Access:**
- Measures freely available for non-commercial use
- Item parameters published
- Administration via REDCap, Epic, or HealthMeasures Assessment Center

**Scraping Feasibility:** HIGH (Legitimate Access)
- Measures downloadable
- Scoring algorithms available
- No scraping needed - official access

**Legal Considerations:** LOW RISK
- NIH-funded public resource
- Free for non-commercial research
- Attribution required

**References:**
- [PROMIS HealthMeasures](https://www.healthmeasures.net/explore-measurement-systems/promis)
- [NIH PROMIS Overview](https://commonfund.nih.gov/promis)

---

### 4.2 PatientsLikeMe

**URL:** https://www.patientslikeme.com/

**Content Type:**
- Self-reported patient health data
- Treatment efficacy reports
- Symptom tracking
- Side effect reports
- Quality of life measures

**Coverage:**
- 830,000+ members
- 2,900+ conditions
- Major focus: ALS, MS, epilepsy, Parkinson's

**Data Structure:**
```
Patient Data:
- Demographics (age, gender, race, ethnicity)
- Conditions with severity
- Symptoms with tracking
- Treatments with efficacy ratings
- Side effects
- Hospitalizations
- Lab results (optional)
```

**Research Access:**
- Partnership/licensing model
- Owned by UnitedHealth Group (Optum Ventures) since 2020
- Research FAQ: https://www.patientslikeme.com/research/faq
- "Data for Good" program available

**Scraping Feasibility:** LOW
- Proprietary data
- Behind authentication
- Corporate data policies

**Legal Considerations:** HIGH RISK
- HIPAA-relevant patient data
- Corporate ownership (UnitedHealth)
- Licensing required

**Data Limitations:**
- Skews female, younger, white, higher education
- Variable activity by condition
- Few seniors (75+)

**Recommendation:** Explore research partnership via Data for Good program. Do not attempt scraping.

**References:**
- [PatientsLikeMe Research](https://www.patientslikeme.com/research/dataforgood)
- [Research FAQs](https://www.patientslikeme.com/research/faq)

---

### 4.3 Open Humans (PRIMARY - Open Access)

**URL:** https://www.openhumans.org/

**Content Type:**
- Participant-donated health data
- Genetic data (23andMe, AncestryDNA, VCF)
- Wearable data (Fitbit, Apple Health)
- Self-reported health information

**Supported Data Sources:**
- 23andMe
- AncestryDNA
- Fitbit
- Runkeeper
- Withings
- uBiome (archived)
- Generic VCF importer
- Apple HealthKit (community bridge)
- FamilyTreeDNA
- Gencove
- Twitter
- Nightscout (diabetes CGM)

**Data Structure:**
- Participant-controlled data files
- Standardized where possible
- Research project-specific formatting

**Access Model:**
- Participant-centered consent
- Researchers create "projects"
- Participants opt-in to share
- API available for approved projects

**Active Projects:**
- Quantified Flu (wearable illness prediction)
- Quality of Life Technologies Lab studies
- Various academic research projects

**Scraping Feasibility:** HIGH (Legitimate)
- Open-source platform
- Research project creation available
- API for approved researchers

**Legal Considerations:** LOW RISK
- 501(c)(3) nonprofit
- Explicit participant consent model
- Open science mission
- Funded by Robert Wood Johnson, Knight, Shuttleworth Foundations

**Recommendation:** HIGHLY RECOMMENDED. Create research project for legitimate data access.

**References:**
- [Open Humans About](https://www.openhumans.org/about/)
- [Open Humans Research](https://www.openhumans.org/create/)
- [GigaScience Paper](https://academic.oup.com/gigascience/article/8/6/giz076/5523201)

---

## 5. Health App Data Sources

### 5.1 Apple HealthKit

**Documentation:** https://developer.apple.com/healthkit/

**Content Type:**
- Steps, distance, activity
- Heart rate, HRV
- Sleep analysis
- Workouts
- Blood pressure, glucose
- Respiratory rate
- Menstrual tracking

**Data Structure:**
- HKQuantitySample (numeric measurements)
- HKCategorySample (discrete values)
- HKWorkout (exercise sessions)
- All timestamps in local time

**Access Model:**
- **NO PUBLIC CLOUD API**
- Local device storage only
- App-by-app permission grants
- Encrypted when device locked
- Requires iOS app development

**Integration Options:**
1. Build native iOS app with HealthKit entitlement
2. Use ResearchKit for research studies
3. Use unified APIs (Terra, Open Wearables) via user sync

**Scraping Feasibility:** LOW (Direct)
- No web scraping possible
- Requires iOS app
- User-by-user consent

**Legal Considerations:** MEDIUM
- Apple Developer Program membership required
- Privacy Review for HealthKit apps
- App Store review guidelines
- HIPAA considerations for health apps

**Reference:** [Apple HealthKit Documentation](https://developer.apple.com/documentation/healthkit)

---

### 5.2 Google Health Connect (Replaces Google Fit)

**Documentation:** https://developer.android.com/health-and-fitness/guides/health-connect

**IMPORTANT: Google Fit Deprecated**
- Google Fit APIs deprecated in 2026
- No new signups since May 1, 2024
- Google Fit officially deprecated July 1, 2025
- Migration to Health Connect required

**Content Type:**
- Activity, exercise, nutrition
- Sleep, vitals, body measurements
- Cycle tracking
- On-device data storage

**Health Connect Features:**
- More secure on-device storage
- Broader device compatibility
- Unified Android health data layer
- Granular permission controls

**Data Structure:**
- Similar data types to Google Fit
- On-device SQLite database
- Permission-based access

**Access Model:**
- Android app required
- User grants specific permissions
- No cloud API (unlike old Google Fit REST API)

**Scraping Feasibility:** LOW (Direct)
- On-device only
- No web interface
- Requires Android app

**Legal Considerations:** MEDIUM
- Google Play Developer policies
- Health data sensitivity requirements

**Reference:** [Google Fit Deprecation Notice](https://developers.google.com/fit)

---

### 5.3 Unified Health APIs (RECOMMENDED)

#### 5.3.1 Terra API

**URL:** https://tryterra.co/

**Supported Devices:**
- Garmin, Fitbit, Apple, Google
- Polar, Eight Sleep, Oura, Whoop
- 100+ devices total

**Features:**
- Single API for all wearables
- OAuth handling
- Data normalization
- Real-time webhooks

**Pricing:** Commercial (contact for research pricing)

#### 5.3.2 Open Wearables (Open Source)

**URL:** https://www.openwearables.io/

**Supported Devices:**
- Garmin, Fitbit, Oura, Whoop
- Apple Health, Strava
- 200+ devices

**Features:**
- Self-hosted platform
- Single API
- OAuth, normalization, syncing
- Open-source

**Integration Time:**
- Days (vs weeks per device manually)
- Handles OAuth, normalization automatically

**Scraping Feasibility:** HIGH
- Open-source self-hosted option
- Commercial unified APIs available

**Legal Considerations:** LOW-MEDIUM
- Aggregator platform, not direct scraping
- User consent required
- Per-device ToS compliance

**References:**
- [Open Wearables](https://www.openwearables.io/)
- [Terra API](https://tryterra.co/)
- [Open Wearables Medium Article](https://medium.com/@themomentum_ai/open-wearables-real-world-use-cases-for-wearable-data-integration-735a85279501)

---

## 6. Wearable / Quantified Self Data Sources

### 6.1 Fitabase (Research Platform - PRIMARY)

**URL:** https://fitabase.com/

**Content Type:**
- Fitbit and Garmin data
- Research-grade data export
- Participant management

**Data Export:**
- CSV files
- Daily, hourly, minute-level data
- Individual or batch export

**Track Record:**
- 1,100+ research studies worldwide
- Academic research focus

**Features:**
- Study management dashboard
- Participant compliance monitoring
- Automated data collection

**Scraping Feasibility:** HIGH (Legitimate Access)
- Designed for researchers
- API and export tools
- Commercial research platform

**Legal Considerations:** LOW RISK
- Research-focused platform
- Proper consent workflows
- Participant agreements

**Pricing:** Commercial (research pricing available)

**Reference:** [Fitabase](https://fitabase.com/)

---

### 6.2 Device-Specific APIs

#### Fitbit Web API
- OAuth 2.0 authentication
- Rate limits apply
- Team approval needed for full data access
- Data types: Activity, sleep, heart rate, SpO2

#### Garmin Connect API
- OAuth 1.0 (less secure)
- Application process required
- Server-to-server integration
- Complex data formats

#### Oura Cloud API
- Sleep, recovery, readiness metrics
- HRV, resting HR, temperature
- Team approval for full access
- Good accuracy reputation

#### Whoop API
- Strain, recovery, sleep
- Developer program required
- Limited public documentation

**Scraping Feasibility:** MEDIUM
- APIs available but access restricted
- Approval processes required
- Rate limits and quotas

**Legal Considerations:** MEDIUM
- Per-platform Terms of Service
- Developer agreement compliance
- User consent requirements

---

### 6.3 Wearipedia Project

**URL:** https://www.medrxiv.org/content/10.1101/2025.05.12.25327465v1

**Content Type:**
- Wearable device documentation
- Clinical research usage
- Privacy/security evaluations
- Open-source coding tools

**Purpose:**
- Help researchers select optimal wearables
- Document data extraction methods
- Track clinical trial usage

**Reference:** [Wearipedia Preprint](https://www.medrxiv.org/content/10.1101/2025.05.12.25327465v1.full.pdf)

---

## 7. Additional Data Sources

### 7.1 WHO GHO OData API

**URL:** https://www.who.int/data/gho/info/gho-odata-api

**Content Type:**
- Global health statistics
- Country-level health indicators
- Disease burden data

**Example Query:**
```
https://ghoapi.azureedge.net/api/DIMENSION/COUNTRY/DimensionValues
```

**Scraping Feasibility:** HIGH
- Open OData API
- No authentication required
- Well-documented

**Legal Considerations:** LOW RISK
- WHO public data
- Attribution required

---

### 7.2 Our World in Data

**URL:** https://ourworldindata.org/

**GitHub:** https://github.com/owid

**Content Type:**
- Global health visualizations
- COVID-19 data
- Health statistics from UN, World Bank

**Access:**
- Charts embeddable with attribution
- Data downloadable
- R package: `owidapi`

**Reference:** [Global Health Data Explorer](https://ourworldindata.org/explorers/global-health)

---

## 8. Scraping Feasibility Summary

### Recommended for Scraping/API Access

| Source | Method | Effort | Notes |
|--------|--------|--------|-------|
| ClinicalTrials.gov | REST API | Low | Primary clinical trials source |
| AACT Database | PostgreSQL dump | Low | Best for bulk analysis |
| WHO ICTRP | Request form | Medium | Global coverage |
| PROMIS | Official download | Low | Free instruments |
| Open Humans | Research project | Medium | Consent-based |
| Open Wearables | Self-host | Medium | Open-source |
| Fitabase | Commercial API | Low | Research-focused |

### Not Recommended (High Legal Risk)

| Source | Reason | Alternative |
|--------|--------|-------------|
| Reddit | Active litigation, ToS | Official Research API application |
| ConsumerLab | Paywall, commercial | NIH ODS API |
| Labdoor | Proprietary scores | PubMed literature |
| Examine.com | Subscription, copyright | NIH ODS, PubMed |
| PatientsLikeMe | Corporate, HIPAA | Open Humans, PROMIS |
| HealthUnlocked | Privacy, ToS | Open Humans |

---

## 9. Legal Framework Summary

### CFAA (Computer Fraud and Abuse Act) Considerations
- Unauthorized access to computer systems is illegal
- "Authorization" increasingly includes ToS violations
- hiQ v. LinkedIn: Public data may be scrapeable, but evolving

### DMCA Anti-Circumvention (Section 1201)
- Reddit lawsuit (Oct 2025) uses this against scrapers
- Bypassing access controls is separate violation
- Applies even to "public" content with technical barriers

### HIPAA Considerations
- Patient health information requires special handling
- De-identification requirements
- Business Associate Agreements may be needed

### GDPR (EU Data)
- Right to data portability
- Consent requirements
- Cross-border transfer restrictions

### Best Practices
1. Always check Terms of Service
2. Use official APIs when available
3. Respect robots.txt
4. Implement rate limiting
5. Store data securely
6. Anonymize before analysis
7. Get legal review for commercial use

---

## 10. Implementation Priority

### Phase 1 - Low Risk, High Value (Immediate)
1. **ClinicalTrials.gov API** - Clinical trial protocols and results
2. **AACT Database** - Bulk clinical trials analysis
3. **PROMIS** - Standardized outcome measures
4. **Open Humans** - Participant-donated multi-source data
5. **WHO GHO API** - Global health statistics

### Phase 2 - Medium Effort (Month 2-3)
1. **WHO ICTRP** - Apply for bulk download access
2. **EU CTIS** - Monitor for API development
3. **Open Wearables** - Self-host for wearable data
4. **Fitabase** - Evaluate for research studies

### Phase 3 - Partnership Required (Month 4+)
1. **Reddit** - Apply via r/reddit4researchers
2. **PatientsLikeMe** - Explore Data for Good program
3. **Commercial APIs** - Terra, device-specific as needed

---

## References

### Clinical Trials
- [ClinicalTrials.gov API](https://clinicaltrials.gov/data-api/api)
- [AACT Database](https://aact.ctti-clinicaltrials.org/)
- [WHO ICTRP](https://www.who.int/tools/clinical-trials-registry-platform)
- [EU CTIS](https://euclinicaltrials.eu/)

### Patient Outcomes
- [PROMIS HealthMeasures](https://www.healthmeasures.net/)
- [Open Humans](https://www.openhumans.org/)
- [PatientsLikeMe Research](https://www.patientslikeme.com/research/dataforgood)

### Health Apps & Wearables
- [Open Wearables](https://www.openwearables.io/)
- [Terra API](https://tryterra.co/)
- [Fitabase](https://fitabase.com/)
- [Apple HealthKit](https://developer.apple.com/healthkit/)
- [Health Connect](https://developer.android.com/health-and-fitness/guides/health-connect)

### Forums & Communities
- [Reddit Developer Platform](https://support.reddithelp.com/hc/en-us/articles/14945211791892-Developer-Platform-Accessing-Reddit-Data)
- [HealthUnlocked](https://healthunlocked.com/)

### Statistics
- [WHO GHO OData API](https://www.who.int/data/gho/info/gho-odata-api)
- [Our World in Data](https://ourworldindata.org/)

---

*Document generated: January 2026*
*For Gene Platform Data Integration*
