# Data Access: Legal Considerations & Governance

**Last Updated:** January 2026
**Purpose:** Actionable guidance for legitimate data access and legal risk mitigation

---

## Executive Summary

| Risk Level | Approach | Examples |
|------------|----------|----------|
| **LOW** | Direct access via APIs or partnerships | Clinical Trials.gov, AACT, PROMIS, Open Humans |
| **MEDIUM** | Restricted access with consent/licensing | Wearable APIs, Health Connect, Fitabase |
| **HIGH** | Avoid or partnership-only | Reddit, ConsumerLab, Labdoor, HealthUnlocked |

**Key Principle:** Use official APIs, respect Terms of Service, and obtain proper consent.

---

## 1. LOW-RISK LEGITIMATE ACCESS

### 1.1 Clinical Trial Data (PRIMARY)

**ClinicalTrials.gov REST API**
- Status: Public domain, NIH-operated
- Access: No registration required
- Data: Protocols, results, enrollment, interventions
- Rate: Unlimited, well-documented
- Compliance: Must follow Terms and Conditions
- References: https://clinicaltrials.gov/data-api/api

**AACT Database** (Recommended for bulk analysis)
- Status: Free derivative of ClinicalTrials.gov
- Access: Free registration, cloud PostgreSQL access
- Data: 50+ normalized tables, daily updates
- Rate: 10 concurrent connections per account
- Compliance: Public domain, cite source
- References: https://aact.ctti-clinicaltrials.org/

**WHO ICTRP** (International coverage)
- Status: WHO public health mandate
- Access: Search portal + bulk download via form request
- Data: Global trials from 17+ registries
- Terms: Attribute source, non-commercial use, keep data current
- References: https://trialsearch.who.int

### 1.2 Patient-Reported Outcomes

**PROMIS** (Standardized measures)
- Status: NIH-funded public resource
- Access: Measures freely downloadable
- Data: Item banks, scoring algorithms, calibrations
- Compliance: Non-commercial research, credit NIH
- Coverage: Pain, fatigue, emotional distress, physical function, sleep
- References: https://www.healthmeasures.net/explore-measurement-systems/promis

**Open Humans** (Participant-consented data)
- Status: 501(c)(3) nonprofit, open science mission
- Access: Create research project, participants opt-in
- Data: Genetic (23andMe, AncestryDNA), wearable (Fitbit, Apple), self-reported
- Compliance: Explicit participant consent, open-source platform
- Advantage: Data donation model aligns with privacy expectations
- References: https://www.openhumans.org/

### 1.3 Global Health Statistics

**WHO GHO OData API**
- Status: WHO public data, open OData API
- Access: No authentication required
- Data: Country-level health indicators, disease burden
- Compliance: Attribution required
- References: https://www.who.int/data/gho/info/gho-odata-api

**Our World in Data**
- Status: Curated global health data
- Access: GitHub repository, downloadable datasets
- Data: COVID-19, health statistics from UN/World Bank
- Compliance: Charts embeddable with attribution
- References: https://ourworldindata.org/

---

## 2. MEDIUM-RISK: CONSENT & LICENSING REQUIRED

### 2.1 Wearable Data

**Fitabase** (Research platform)
- Status: Commercial, research-focused
- Access: Designed for researchers with proper consent
- Data: Fitbit, Garmin research-grade exports
- Process: Study management, participant compliance monitoring
- Compliance: Proper consent workflows, participant agreements
- Track Record: 1,100+ research studies worldwide

**Open Wearables** (Self-hosted)
- Status: Open-source unified API
- Access: Self-hosted or commercial platform
- Devices: Garmin, Fitbit, Oura, Whoop, Apple, Strava (200+)
- Compliance: Per-device ToS compliance, user consent required
- Advantage: Single API normalizes across devices

**Device-Specific APIs** (Fitbit, Garmin, Oura, Whoop)
- Status: OAuth 2.0 authentication (most)
- Access: Developer account + approval process required
- Compliance: Per-platform Terms of Service, rate limits, quotas
- Considerations: User-by-user consent, server-to-server integration

### 2.2 Native Mobile Health Data

**Apple HealthKit**
- Status: Local device storage only, no public cloud API
- Access: Requires iOS app development with HealthKit entitlement
- Compliance: Apple Developer Program membership, Privacy Review, App Store guidelines
- Data Types: Steps, heart rate, sleep, blood pressure, glucose, respiratory rate
- Note: User-by-user consent, encrypted when locked

**Google Health Connect** (Replaces deprecated Google Fit)
- Status: On-device storage, no cloud API
- Access: Android app required with permission grants
- Compliance: Google Play Developer policies
- Data Types: Activity, sleep, vitals, body measurements, cycle tracking
- Timeline: Google Fit deprecated July 2025; must migrate to Health Connect

---

## 3. HIGH-RISK: AVOID OR PARTNERSHIP ONLY

### 3.1 Why These Are High Risk

**Legal Exposure:**
- Terms of Service violations can trigger CFAA (Computer Fraud and Abuse Act) liability
- DMCA Section 1201 anti-circumvention claims (Reddit lawsuit vs. Perplexity, October 2025)
- Proprietary data protection and copyright claims
- HIPAA violations for unprotected health data

**Problematic Sources:**

| Source | Risk | Reason | Alternative |
|--------|------|--------|-------------|
| **Reddit** | HIGH | Active litigation, ToS prohibits scraping, DMCA claims | Apply via r/reddit4researchers for official Research API |
| **HealthUnlocked** | HIGH | Medical data privacy, patient expectation of privacy, ToS restrictions | Partnership approach or Open Humans |
| **ConsumerLab** | HIGH | Paywall, commercial entity, copyrighted test results | NIH ODS API for supplement science |
| **Labdoor** | HIGH | Proprietary scoring, VC-backed, commercial entity | PubMed literature for supplement evidence |
| **Examine.com** | HIGH | Subscription model, copyrighted editorial content | NIH ODS, PubMed for evidence-based data |
| **PatientsLikeMe** | HIGH | Corporate ownership (UnitedHealth), HIPAA-relevant data, licensing required | Explore "Data for Good" partnership program |

### 3.2 Reddit Case Study: Why Official Access Matters

**What Happened:**
- Reddit sued Perplexity AI (October 2025) for DMCA circumvention
- Named third-party scrapers: SerpApi, Oxylabs
- Claims: Unauthorized access despite public content

**Correct Approach:**
1. Apply via r/reddit4researchers for formal Research API access
2. Submit research proposal with clear methodology
3. Follow Responsible Builder Policy
4. Respect rate limits (60 requests/minute)
5. Anonymize data in publications, credit Reddit

**References:**
- https://support.reddithelp.com/hc/en-us/articles/14945211791892-Developer-Platform-Accessing-Reddit-Data
- https://support.reddithelp.com/hc/en-us/articles/42728983564564-Responsible-Builder-Policy

---

## 4. LEGAL FRAMEWORK SUMMARY

### Computer Fraud and Abuse Act (CFAA)
- Unauthorized access to computer systems is illegal
- "Authorization" increasingly includes ToS compliance
- Scraping despite explicit ToS prohibitions may constitute violation
- Individual counts per violation can compound liability

### DMCA Section 1201 (Anti-Circumvention)
- Bypassing access controls, even for public data, is separate violation
- Applies independently of copyright infringement
- Recently used against web scrapers (Reddit v. Perplexity, 2025)
- Circumvention includes: account takeover, credential re-use, ToS evasion

### HIPAA (Health Information Privacy)
- Applies to Protected Health Information (PHI) in identifiable form
- De-identification standards exist (HIPAA Safe Harbor, Expert Determination)
- Business Associate Agreements may be required for data processing
- Most important for: patient health records, condition-sensitive data

### GDPR (EU Data Protection)
- Right to data portability (but only for consenting individuals)
- Explicit consent required for processing
- Cross-border transfer restrictions (EU → non-EU)
- Important if analyzing EU residents' health data

### Copyright & Database Rights
- Published research summaries (Examine.com, ConsumerLab) are copyrighted
- Database designs themselves may have legal protection
- Fair use may not apply to commercial product alternatives

---

## 5. DECISION FRAMEWORK

### Before Accessing Any Data Source:

**Step 1: Check Official Access**
- [ ] Does the source offer an official API?
- [ ] Is there a free tier or research license?
- [ ] Can I register for legitimate access?
- **Outcome:** Use official access if available

**Step 2: Review Terms of Service**
- [ ] Does ToS explicitly prohibit scraping?
- [ ] Are there researcher exceptions?
- [ ] What are rate limits and quotas?
- [ ] Is there a research approval process?
- **Outcome:** Determine compliance path

**Step 3: Assess Data Sensitivity**
- [ ] Is this health/medical data (HIPAA)?
- [ ] Contains personally identifiable information (PII)?
- [ ] Includes genetic data?
- [ ] Are EU residents' data included (GDPR)?
- **Outcome:** Identify compliance requirements

**Step 4: Evaluate Technical Barriers**
- [ ] Content is publicly visible without login?
- [ ] Or requires authentication/account?
- [ ] Or requires technical circumvention?
- **Outcome:** Determine CFAA/DMCA risk

**Step 5: Make Decision**
- **GREEN:** Official API + no sensitivity issues → Proceed
- **YELLOW:** Restricted access available → Pursue partnership/licensing
- **RED:** High sensitivity + no official access → Do not attempt scraping

---

## 6. RECOMMENDED IMPLEMENTATION PRIORITY

### Phase 1: Low Risk, Immediate Implementation
1. **ClinicalTrials.gov API** - Clinical trial protocols and results
2. **AACT Database** - Bulk clinical trial analysis
3. **PROMIS** - Standardized patient-reported outcome measures
4. **Open Humans** - Participant-donated multi-source genetic/wearable data
5. **WHO GHO API** - Global health statistics

**Expected Timeline:** 2-4 weeks
**Legal Review Needed:** Minimal (all public domain/open science)

### Phase 2: Medium Effort (Month 2-3)
1. **WHO ICTRP** - Apply for bulk download access (form-based)
2. **EU CTIS** - Monitor API development
3. **Open Wearables** - Self-host for wearable data aggregation
4. **Fitabase** - Evaluate research licensing for study designs

**Expected Timeline:** 4-8 weeks
**Legal Review Needed:** Standard (review ToS, confirm compliance)

### Phase 3: Partnership-Based (Month 4+)
1. **Reddit** - Apply via r/reddit4researchers official process
2. **PatientsLikeMe** - Explore "Data for Good" partnership program
3. **Commercial Wearable APIs** - Negotiate research rates (Terra, device-specific)

**Expected Timeline:** 8-16 weeks
**Legal Review Needed:** Full (partnership agreements, data processing)

---

## 7. BEST PRACTICES

### Before Launch
- [ ] Legal review of all data sources in use
- [ ] Document why each source was chosen
- [ ] Maintain approval records (API keys, partnership agreements)
- [ ] Review HIPAA compliance if health data involved

### During Data Collection
- [ ] Implement rate limiting (respect API quotas)
- [ ] Log all data access (audit trail)
- [ ] Store data securely (encrypted at rest)
- [ ] Document data lineage (source, date, version)

### Before Analysis/Publication
- [ ] Anonymize sensitive data (remove PII, genetic identifiers)
- [ ] Credit data sources in publications
- [ ] Cite original researchers and databases
- [ ] Confirm ToS allow publication of results

### Continuous Compliance
- [ ] Monitor source ToS changes (set calendar reminders)
- [ ] Check for API deprecations (e.g., Google Fit → Health Connect)
- [ ] Update partnership agreements as needed
- [ ] Track regulatory changes (HIPAA, GDPR, CFAA interpretation)

---

## 8. ESCALATION CONTACTS

**When to Seek Legal Review:**
- Considering scraping any source not explicitly approved
- Processing health data from EU residents (GDPR)
- Involving patient data with identifiers (HIPAA)
- Commercial product development using third-party data
- Disagreement about Terms of Service interpretation

**Recommended:** Consult data governance legal counsel before:
- Launching new data pipeline
- Processing sensitive health/genetic data
- International data transfers
- Commercial product launches

---

## References

### Clinical Trials
- [ClinicalTrials.gov API](https://clinicaltrials.gov/data-api/api)
- [AACT Database](https://aact.ctti-clinicaltrials.org/)
- [WHO ICTRP](https://www.who.int/tools/clinical-trials-registry-platform)

### Patient Outcomes
- [PROMIS HealthMeasures](https://www.healthmeasures.net/)
- [Open Humans](https://www.openhumans.org/)
- [PatientsLikeMe Data for Good](https://www.patientslikeme.com/research/dataforgood)

### Health Apps & Wearables
- [Open Wearables](https://www.openwearables.io/)
- [Terra API](https://tryterra.co/)
- [Fitabase](https://fitabase.com/)
- [Apple HealthKit](https://developer.apple.com/healthkit/)
- [Health Connect](https://developer.android.com/health-and-fitness/guides/health-connect)

### Legal References
- CFAA: 18 U.S.C. § 1030
- DMCA: 17 U.S.C. § 1201
- HIPAA: 45 CFR Parts 160 and 164
- GDPR: Regulation 2016/679

---

**Next Steps:** Reference this document when planning data integration. For specific sources, confirm with legal team before implementation.
