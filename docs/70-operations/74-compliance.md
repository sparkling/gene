# Regulatory Compliance

**Document ID:** 74-COMPLIANCE
**Status:** Final
**Owner:** Legal/Compliance
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Platform positioned as educational/informational tool, NOT medical device. Key compliance requirements: (1) FDA - avoid SaMD classification through educational positioning and disclaimers, (2) FTC - substantiated health claims only, (3) HIPAA - privacy-first architecture with field-level encryption for genetic data, (4) GDPR - explicit consent, data portability, right to deletion. Target: SOC 2 Type II certification. Medical disclaimers required on all recommendations.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| FDA positioning | Educational, not diagnostic | Avoid SaMD regulation | Jan 2026 |
| Health claims approach | Evidence-linked only | FTC compliance | Jan 2026 |
| Data encryption | Field-level AES-256 | HIPAA-ready | Jan 2026 |
| Consent model | Explicit opt-in | GDPR compliance | Jan 2026 |
| Compliance target | SOC 2 Type II | Enterprise readiness | Jan 2026 |

---

## Regulatory Landscape

### Applicable Regulations

| Regulation | Jurisdiction | Applicability | Priority |
|------------|--------------|---------------|----------|
| FDA 21 CFR Part 820 | US | Software as Medical Device | HIGH |
| FTC Act Section 5 | US | Health claims, advertising | HIGH |
| HIPAA | US | Health data privacy | HIGH |
| GINA | US | Genetic non-discrimination | HIGH |
| GDPR | EU | Data privacy | HIGH |
| CCPA/CPRA | California | Data privacy | MEDIUM |
| PIPEDA | Canada | Data privacy | MEDIUM |
| LGPD | Brazil | Data privacy | LOW |

---

## FDA Software as Medical Device

### Classification Analysis

#### What is a Medical Device?

FDA defines a medical device as an instrument intended for use in:
- Diagnosis of disease
- Cure, mitigation, treatment, or prevention of disease
- Affecting structure or function of the body

#### SaMD (Software as Medical Device) Risk Categories

| Category | Risk Level | Intended Use | Our Status |
|----------|------------|--------------|------------|
| Class I | Low | General wellness, education | **TARGET** |
| Class II | Moderate | Diagnosis support | AVOID |
| Class III | High | Treatment decisions | AVOID |

### Our Positioning: General Wellness / Education

**Key Distinction:** We provide INFORMATION, not DIAGNOSIS or TREATMENT RECOMMENDATIONS.

```
WHAT WE ARE:                        WHAT WE ARE NOT:
├── Educational platform            ├── Diagnostic tool
├── Information aggregator          ├── Treatment prescriber
├── Research companion              ├── Medical device
├── Knowledge connector             ├── Clinical decision support
└── Community platform              └── Health intervention system
```

### FDA Compliance Strategy

#### 1. Intended Use Statement

> "[Platform Name] is an educational platform that helps users understand their genetic data in the context of publicly available scientific research. The platform does not diagnose, treat, cure, or prevent any disease. All information is for educational purposes only and should not be considered medical advice. Users should consult qualified healthcare providers before making any health decisions."

#### 2. Feature Design Principles

| Principle | Implementation |
|-----------|----------------|
| No diagnosis language | "Your data shows..." not "You have..." |
| No treatment commands | "Research suggests..." not "You should take..." |
| Evidence attribution | All claims linked to published research |
| User agency | "Explore" and "learn" not "treat" |
| Professional referral | Encourage consulting healthcare providers |

#### 3. Content Guidelines

**Allowed:**
- "Research has associated this variant with..."
- "Studies suggest this pathway may be involved in..."
- "Some practitioners recommend exploring..."
- "You may want to discuss with your doctor..."

**Not Allowed:**
- "This variant causes your symptoms"
- "Take this supplement to fix..."
- "This will cure your condition"
- "Stop taking your medication because..."

### FDA Exemptions We Rely On

1. **General Wellness Exemption**: Products intended for general wellness (e.g., weight management, physical fitness, relaxation) that present low risk
2. **Educational Exemption**: Products that provide educational information without clinical claims
3. **21st Century Cures Act**: Clinical decision support software exemptions for certain types of software

---

## FTC Health Claims

### FTC Requirements

The FTC requires that:
1. Health claims be **truthful and not misleading**
2. Claims be **substantiated** before being made
3. **Material connections** be disclosed (affiliates, sponsorships)

### Claim Substantiation Framework

| Claim Type | Evidence Required | Our Approach |
|------------|-------------------|--------------|
| Structure/function | Competent scientific evidence | Link to peer-reviewed research |
| Health benefit | RCTs or strong observational | Only cite published studies |
| Cure/treatment | FDA approval | NEVER make these claims |

### Implementation

#### 1. Evidence Grading System

| Grade | Evidence Type | Display |
|-------|---------------|---------|
| A | Multiple RCTs | "Strong evidence" |
| B | Single RCT or meta-analysis | "Moderate evidence" |
| C | Observational studies | "Preliminary evidence" |
| D | Animal/cell studies | "Early research" |
| E | Traditional use only | "Traditional use" |

#### 2. Claim Review Process

```
CONTENT CLAIM REVIEW WORKFLOW
        │
        ▼
┌───────────────┐
│ Writer drafts │
│ content/claim │
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ Evidence check│
│ (researcher)  │
└───────┬───────┘
        │
        ▼
┌───────────────┐     ┌───────────────┐
│ Legal review  │────▶│ Compliance    │
│ (if flagged)  │     │ database      │
└───────┬───────┘     └───────────────┘
        │
        ▼
┌───────────────┐
│ Publication   │
└───────────────┘
```

#### 3. Affiliate Disclosure

For any affiliate relationships (supplements, lab tests):
- Clear disclosure near the recommendation
- "We may earn a commission" language
- No influence on recommendations

---

## HIPAA Requirements

### Applicability

HIPAA applies to:
- **Covered Entities**: Healthcare providers, health plans, clearinghouses
- **Business Associates**: Those handling PHI on behalf of covered entities

**Our Status:** We are NOT a covered entity, but we handle health data. We adopt HIPAA-grade security as best practice and to enable future B2B relationships.

### Technical Safeguards

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| Access controls | Role-based access (RBAC) | Planned |
| Audit controls | Comprehensive logging | Planned |
| Integrity controls | Hash verification | Planned |
| Transmission security | TLS 1.3 | Implemented |
| Encryption at rest | AES-256 | Implemented |

### Administrative Safeguards

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| Security officer | Designated role | Planned |
| Workforce training | Annual HIPAA training | Planned |
| Access management | Principle of least privilege | Planned |
| Incident response | Documented procedures | Planned |
| Risk assessment | Annual assessment | Planned |

### Physical Safeguards

| Requirement | Implementation | Status |
|-------------|----------------|--------|
| Facility access | Cloud provider controls (AWS/GCP) | Inherited |
| Workstation use | Remote work policy | Planned |
| Device controls | Encryption required | Planned |

### Data Handling

#### Genetic Data Classification

| Data Type | Sensitivity | Protection |
|-----------|-------------|------------|
| Raw genetic file | HIGH | Field-level encryption |
| SNP interpretations | HIGH | Field-level encryption |
| Lab results | HIGH | Field-level encryption |
| Symptoms | MEDIUM | Database encryption |
| Profile info | MEDIUM | Database encryption |
| Usage analytics | LOW | Standard protection |

#### Encryption Architecture

```
USER DATA FLOW
        │
        ▼
┌───────────────┐
│ Client-side   │ ◄── TLS 1.3 in transit
│ (Browser)     │
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ API Gateway   │ ◄── JWT authentication
│               │
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ Application   │ ◄── Field-level encryption before storage
│ Layer         │
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ PostgreSQL    │ ◄── AES-256 at rest
│ (encrypted)   │ ◄── Encrypted backups
└───────────────┘
```

---

## GDPR Requirements

### Applicability

GDPR applies when:
- Processing data of EU residents
- Offering services to EU market
- Monitoring behavior of EU individuals

**Our Status:** Processing genetic data (special category) requires explicit consent.

### Lawful Basis

| Basis | Applicability | Our Use |
|-------|---------------|---------|
| Consent | Special category data | PRIMARY |
| Contract | Service delivery | SECONDARY |
| Legitimate interest | Analytics | LIMITED |

### User Rights Implementation

| Right | Implementation | Timeline |
|-------|----------------|----------|
| Right to access | Data export feature | Immediate |
| Right to rectification | Profile editing | Immediate |
| Right to erasure | Account deletion | 30 days |
| Right to portability | JSON/CSV export | Immediate |
| Right to restrict | Processing pause | Immediate |
| Right to object | Marketing opt-out | Immediate |

### Consent Requirements

#### Consent for Genetic Data

```
CONSENT FLOW
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ GENETIC DATA CONSENT                                   │
│                                                        │
│ We will process your genetic data to:                  │
│ □ Analyze your uploaded DNA file                       │
│ □ Provide personalized health insights                 │
│ □ Connect your data to our knowledge database          │
│                                                        │
│ This data is:                                          │
│ • Encrypted with AES-256 encryption                   │
│ • Never sold to third parties                         │
│ • Deletable at any time via account settings          │
│                                                        │
│ [I CONSENT - Required to use service]                  │
│ [LEARN MORE ABOUT DATA HANDLING]                       │
└───────────────────────────────────────────────────────┘
```

#### Consent for Research (Optional)

```
OPTIONAL: CONTRIBUTE TO RESEARCH
        │
        ▼
┌───────────────────────────────────────────────────────┐
│ ANONYMOUS RESEARCH CONTRIBUTION                        │
│                                                        │
│ Help advance scientific understanding by allowing      │
│ your anonymized, de-identified data to be included    │
│ in aggregate research datasets.                        │
│                                                        │
│ • Your identity is NEVER included                     │
│ • Data is irreversibly anonymized                     │
│ • You can withdraw consent at any time                │
│                                                        │
│ □ Yes, I want to contribute to research               │
│ □ No, keep my data private                            │
└───────────────────────────────────────────────────────┘
```

### Data Processing Records

Maintain records of:
- Categories of data processed
- Purposes of processing
- Data retention periods
- Security measures
- Third-party processors

---

## GINA (Genetic Information Nondiscrimination Act)

### Requirements

GINA prohibits:
- Health insurers from using genetic information for coverage/pricing
- Employers from using genetic information for employment decisions

### Our Obligations

1. **Never share genetic data** with insurers or employers
2. **Educate users** about their GINA protections
3. **Document compliance** in privacy policy

### User Education

Include in onboarding:
> "Your genetic data is protected by GINA (Genetic Information Nondiscrimination Act). Health insurers and employers cannot use your genetic information to discriminate against you. We will never share your genetic data with insurers or employers."

---

## SOC 2 Compliance Target

### SOC 2 Trust Service Criteria

| Criterion | Description | Priority |
|-----------|-------------|----------|
| Security | Protection against unauthorized access | HIGH |
| Availability | System availability for operation | HIGH |
| Processing Integrity | Complete, accurate processing | MEDIUM |
| Confidentiality | Protection of confidential data | HIGH |
| Privacy | Personal information protection | HIGH |

### SOC 2 Roadmap

| Phase | Timeline | Activities |
|-------|----------|------------|
| Gap Assessment | Q2 2026 | Identify control gaps |
| Remediation | Q3 2026 | Implement missing controls |
| Type I Audit | Q4 2026 | Point-in-time assessment |
| Monitoring | Q1-Q2 2027 | 6-month observation |
| Type II Audit | Q3 2027 | Period-of-time assessment |

---

## Medical Disclaimer Framework

### Standard Disclaimer (Footer/Terms)

> **Medical Disclaimer**
>
> [Platform Name] is an educational platform providing information about genetics, nutrition, and traditional medicine. The content on this platform is for informational and educational purposes only and is not intended as medical advice, diagnosis, or treatment.
>
> Always seek the advice of your physician or other qualified healthcare provider with any questions you may have regarding a medical condition. Never disregard professional medical advice or delay in seeking it because of something you have read on this platform.
>
> Genetic information is just one factor in health. Results should be considered in the context of your complete health picture. The presence or absence of genetic variants does not guarantee the development or absence of any condition.
>
> If you think you may have a medical emergency, call your doctor, go to the emergency department, or call emergency services immediately.

### Contextual Disclaimers

| Context | Disclaimer |
|---------|------------|
| SNP report | "This information is educational. Discuss with your healthcare provider before making health decisions." |
| Supplement suggestion | "Traditional use is not a substitute for medical treatment. Consult your doctor before starting any supplement." |
| Herb recommendation | "These herbs may interact with medications. Always consult a qualified practitioner." |
| AI response | "This response is generated from research literature and is not medical advice." |

---

## Compliance Monitoring

### Audit Schedule

| Audit Type | Frequency | Scope |
|------------|-----------|-------|
| Internal compliance | Quarterly | All policies |
| Security assessment | Annually | Technical controls |
| Privacy assessment | Annually | Data handling |
| Legal review | Annually | Terms, disclaimers |
| External audit | Per SOC 2 | Full scope |

### Incident Response

```
COMPLIANCE INCIDENT RESPONSE
        │
        ▼
┌───────────────┐
│ Incident      │ ◄── Identify and log
│ Detection     │
└───────┬───────┘
        │
        ▼
┌───────────────┐
│ Assessment    │ ◄── Severity, scope, impact
│               │
└───────┬───────┘
        │
    ┌───┴───┐
    │       │
    ▼       ▼
┌───────┐ ┌───────┐
│Minor  │ │Major  │
│(log)  │ │(alert)│
└───────┘ └───┬───┘
              │
              ▼
        ┌───────────────┐
        │ Notification  │ ◄── Within 72 hours (GDPR)
        │ (if required) │
        └───────┬───────┘
                │
                ▼
        ┌───────────────┐
        │ Remediation   │
        │ & Post-mortem │
        └───────────────┘
```

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [44-ARCHITECTURE](../40-product/44-architecture.md) | Security architecture |
| [90-RISK-OVERVIEW](../90-risk/90-RISK-OVERVIEW.md) | Regulatory risk |
| [92-THREATS](../90-risk/92-THREATS.md) | Compliance threats |

---

## Open Questions

- [ ] Engage healthcare regulatory counsel for FDA strategy validation
- [ ] Determine timing for SOC 2 Type I audit
- [ ] Assess need for HIPAA Business Associate Agreements with B2B customers
- [ ] Evaluate international data transfer mechanisms (EU-US)

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Legal | Complete compliance framework |
