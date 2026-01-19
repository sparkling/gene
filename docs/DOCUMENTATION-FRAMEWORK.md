# Documentation Framework (McKinsey Standard)

**Document ID:** FRAMEWORK
**Status:** Final
**Owner:** TBD
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

This framework defines how all strategic and operational documentation is structured, organized, and maintained. It follows McKinsey consulting standards and the Pyramid Principle. All documents use a consistent template and adhere to MECE principles.

---

## Governing Principles

### 1. The Pyramid Principle (Barbara Minto)

The foundational McKinsey communication framework:

1. **Start with the answer** — Lead with the recommendation/conclusion
2. **Group supporting arguments** — Cluster related ideas logically
3. **Order within groups** — Use logical flow (time, structure, importance)
4. **Ensure MECE coverage** — Mutually Exclusive, Collectively Exhaustive

**Application:**
- Every document begins with a TL;DR (the answer)
- Supporting sections are grouped by theme
- No gaps in coverage, no overlapping content

---

### 2. The Strategic Choices Cascade

Five questions that define strategy (Lafley/Martin framework):

| Question | Maps To |
|----------|---------|
| What is our winning aspiration? | Vision, Mission |
| Where will we play? | Market, Segments, Geography |
| How will we win? | Value Proposition, Positioning |
| What capabilities must we have? | Product, Technology, Data |
| What management systems are required? | Operations, Processes |

---

### 3. MECE Principle

**Mutually Exclusive, Collectively Exhaustive**

- **Mutually Exclusive:** Categories do not overlap
- **Collectively Exhaustive:** Categories cover everything

**Application:**
- Document sections do not duplicate content
- Together, all documents cover the complete strategic picture
- Each piece of information lives in exactly one place

---

### 4. Document Quality Standards

| Standard | Description |
|----------|-------------|
| Evidence-backed | Every claim supported by evidence or explicit assumption |
| No orphan assertions | Unsupported statements flagged as assumptions |
| Clear ownership | Each document has a designated owner |
| Version control | All changes logged with date and author |
| Cross-referenced | Dependencies between documents explicitly stated |

---

## Document Architecture

```
docs/
│
├── 00-INDEX.md                      # Master navigation
├── 01-EXECUTIVE-SUMMARY.md          # Pyramid Principle: answer first
├── DOCUMENTATION-FRAMEWORK.md       # This document
│
├── 10-strategy/
│   ├── 10-STRATEGY.md               # Strategic foundation
│   ├── 11-VISION-MISSION.md         # Vision, Mission, Values
│   ├── 12-VALUE-PROPOSITION.md      # Value prop deep-dive
│   ├── 13-POSITIONING.md            # Market positioning
│   └── 14-BRAND.md                  # Brand identity, naming
│
├── 20-market/
│   ├── 20-MARKET-OVERVIEW.md        # Market analysis summary
│   ├── 21-MARKET-SIZE.md            # TAM/SAM/SOM analysis
│   ├── 22-MARKET-TRENDS.md          # Trends and dynamics
│   ├── 23-COMPETITORS.md            # Competitive landscape
│   ├── 24-COMPETITOR-PROFILES.md    # Deep competitor analysis
│   └── 25-WHITE-SPACE.md            # Opportunity gaps
│
├── 30-customers/
│   ├── 30-CUSTOMER-OVERVIEW.md      # Customer strategy summary
│   ├── 31-SEGMENTS.md               # Market segmentation
│   ├── 32-PERSONAS.md               # Detailed user personas
│   ├── 33-JOBS-TO-BE-DONE.md        # JTBD framework
│   └── 34-USER-JOURNEYS.md          # Journey maps
│
├── 40-product/
│   ├── 40-PRODUCT-OVERVIEW.md       # Product strategy summary
│   ├── 41-FEATURES.md               # Feature inventory
│   ├── 42-ROADMAP.md                # Product roadmap
│   ├── 43-DATA-SOURCES.md           # THREE WORLDS databases
│   ├── 44-ARCHITECTURE.md           # Technical architecture
│   └── 45-DATA-MODEL.md             # Core data model
│
├── 50-business/
│   ├── 50-BUSINESS-MODEL.md         # Business model summary
│   ├── 51-REVENUE-STREAMS.md        # Revenue model detail
│   ├── 52-PRICING.md                # Pricing strategy
│   ├── 53-MARKETPLACE.md            # Marketplace dynamics
│   └── 54-UNIT-ECONOMICS.md         # Unit economics
│
├── 60-gtm/
│   ├── 60-GTM-OVERVIEW.md           # Go-to-market summary
│   ├── 61-LAUNCH-STRATEGY.md        # Launch plan
│   ├── 62-CHANNELS.md               # Channel strategy
│   ├── 63-MARKETING.md              # Marketing approach
│   └── 64-PARTNERSHIPS.md           # Partnership strategy
│
├── 70-operations/
│   ├── 70-OPERATIONS-OVERVIEW.md    # Operations summary
│   ├── 71-TEAM.md                   # Team structure
│   ├── 72-PROCESSES.md              # Key processes
│   ├── 73-TECHNOLOGY.md             # Technology stack
│   └── 74-COMPLIANCE.md             # Regulatory compliance
│
├── 80-financials/
│   ├── 80-FINANCIALS-OVERVIEW.md    # Financial summary
│   ├── 81-PROJECTIONS.md            # Revenue projections
│   ├── 82-COSTS.md                  # Cost structure
│   ├── 83-FUNDING.md                # Funding requirements
│   └── 84-METRICS.md                # KPIs and metrics
│
└── 90-risk/
    ├── 90-RISK-OVERVIEW.md          # Risk summary
    ├── 91-ASSUMPTIONS.md            # Critical assumptions
    ├── 92-THREATS.md                # Competitive/market threats
    └── 93-MITIGATIONS.md            # Risk mitigations
```

---

## Document Template Standard

Every document follows this structure:

```markdown
# [Document Title]

**Document ID:** [XX-NAME]
**Status:** Draft | Review | Final
**Owner:** [Name/Role]
**Last Updated:** [Date]
**Version:** [X.X]

---

## TL;DR

[Maximum 3 sentences. Answer first per Pyramid Principle. What is the key takeaway?]

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| [What was decided] | [What we chose] | [Why] | [When] |

---

## Context

[Background information needed to understand this document. Why does this document exist? What problem does it address?]

---

## Analysis

[MECE breakdown of the topic. Use headers, tables, and visuals to organize.]

### [Section 1]
...

### [Section 2]
...

---

## Recommendations

[Clear, actionable recommendations. Numbered for easy reference.]

1. **[Recommendation 1]** — [Brief explanation]
2. **[Recommendation 2]** — [Brief explanation]

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [XX-NAME] | Depends on / Informs |

---

## Open Questions

- [ ] [Unresolved question 1]
- [ ] [Unresolved question 2]

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | [Date] | [Author] | Initial version |

---

## Appendix

[Supporting details, raw data, references]
```

---

## Document Specifications

### Root Documents

#### 00-INDEX.md
**Purpose:** Master navigation and document status tracker

**Contains:**
- Document inventory with status
- Reading order recommendations by audience
- Quick links to key sections
- Document dependency map

---

#### 01-EXECUTIVE-SUMMARY.md
**Purpose:** Pyramid Principle lead document — answer first

**Contains:**
- One-page strategic summary
- Key decisions made
- Critical numbers (market size, projections)
- Investment thesis (if applicable)
- Recommended actions

**Structure:**
```
1. The Opportunity (2-3 sentences)
2. Our Solution (2-3 sentences)
3. Why Now (2-3 sentences)
4. Why Us (2-3 sentences)
5. Key Metrics (table)
6. Ask (if applicable)
```

---

### Strategy (10-series)

#### 10-STRATEGY.md
**Purpose:** Strategic foundation consolidation

**Contains:**
- Vision statement with rationale
- Mission statement with scoring
- Value proposition with scoring
- Positioning statement with scoring
- Tagline with scoring
- Strategic choices cascade mapping

---

#### 11-VISION-MISSION.md
**Purpose:** Deep-dive on vision and mission

**Contains:**
- Vision statement analysis
- Mission statement analysis
- McKinsey scoring methodology
- Iteration history
- Final rationale

---

#### 12-VALUE-PROPOSITION.md
**Purpose:** Value proposition deep-dive

**Contains:**
- Full value proposition
- McKinsey VP framework
- Element prioritization matrix
- Scoring and rationale
- Platform offerings inventory

---

#### 13-POSITIONING.md
**Purpose:** Market positioning deep-dive

**Contains:**
- Positioning statement (McKinsey format)
- Segment analysis and decision
- Competitive positioning matrix
- Perceptual map
- THREE WORLDS framework

---

#### 14-BRAND.md
**Purpose:** Brand identity documentation

**Contains:**
- Brand name
- Tagline and rationale
- Brand voice guidelines
- Visual identity principles
- Naming criteria and candidates

---

### Market (20-series)

#### 20-MARKET-OVERVIEW.md
**Purpose:** Market analysis summary

**Contains:**
- Market definition
- Key findings summary
- Market attractiveness assessment
- Strategic implications

---

#### 21-MARKET-SIZE.md
**Purpose:** Total Addressable Market analysis

**Contains:**
- TAM calculation methodology
- SAM (Serviceable) definition
- SOM (Obtainable) targets
- Market growth projections
- Data sources and assumptions

---

#### 22-MARKET-TRENDS.md
**Purpose:** Market dynamics and trends

**Contains:**
- Macro trends affecting the market
- Technology trends
- Consumer behavior trends
- Regulatory trends
- Emerging segments

---

#### 23-COMPETITORS.md
**Purpose:** Competitive landscape overview

**Contains:**
- Competitor tier classification
- Feature comparison matrix
- Positioning map
- Competitive gaps
- Threat assessment

---

#### 24-COMPETITOR-PROFILES.md
**Purpose:** Deep competitor analysis

**Contains:**
- Individual competitor profiles
- Strengths/weaknesses analysis
- Pricing analysis
- Feature deep-dives
- Strategic implications

---

#### 25-WHITE-SPACE.md
**Purpose:** Opportunity gap analysis

**Contains:**
- Unmet needs mapping
- Feature gaps by competitor
- Market positioning gaps
- Technology gaps
- Differentiation opportunities

---

### Customers (30-series)

#### 30-CUSTOMER-OVERVIEW.md
**Purpose:** Customer strategy summary

**Contains:**
- Target market definition
- Segmentation approach
- Primary vs secondary segments
- Customer acquisition strategy overview

---

#### 31-SEGMENTS.md
**Purpose:** Market segmentation analysis

**Contains:**
- Segmentation methodology
- Segment profiles
- Segment sizing
- Segment prioritization
- Selection rationale

---

#### 32-PERSONAS.md
**Purpose:** Detailed user personas

**Contains:**
- 3-5 detailed personas
- Demographics, psychographics
- Goals, frustrations, behaviors
- Technology usage
- Decision-making process

**Persona Template:**
```
Name: [Fictional name]
Segment: [Which segment]
Age/Location: [Demographics]
Condition: [Health situation]
Journey Stage: [Awareness → Consideration → Decision]
Quote: "[In their own words]"
Goals: [What they want to achieve]
Frustrations: [Current pain points]
Behaviors: [How they research, decide]
Channels: [Where they spend time]
Objections: [Why they might not buy]
```

---

#### 33-JOBS-TO-BE-DONE.md
**Purpose:** JTBD framework analysis

**Contains:**
- Functional jobs
- Emotional jobs
- Social jobs
- Job mapping to features
- Outcome expectations

---

#### 34-USER-JOURNEYS.md
**Purpose:** Customer journey mapping

**Contains:**
- Awareness stage touchpoints
- Consideration stage touchpoints
- Decision stage touchpoints
- Onboarding journey
- Ongoing engagement journey
- Pain points and opportunities by stage

---

### Product (40-series)

#### 40-PRODUCT-OVERVIEW.md
**Purpose:** Product strategy summary

**Contains:**
- Product vision
- Core value delivery
- Platform architecture overview
- Feature prioritization approach
- Roadmap summary

---

#### 41-FEATURES.md
**Purpose:** Feature inventory and prioritization

**Contains:**
- Complete feature list
- MVP vs future features
- Feature prioritization matrix
- Build vs buy decisions
- Feature dependencies

---

#### 42-ROADMAP.md
**Purpose:** Product roadmap

**Contains:**
- Phase 1 (MVP) scope
- Phase 2 scope
- Phase 3 scope
- Timeline estimates
- Dependencies and risks

---

#### 43-DATA-SOURCES.md
**Purpose:** THREE WORLDS database inventory

**Contains:**
- Complete data source inventory
- License and access details
- Integration priority
- Data quality assessment
- Update frequency

---

#### 44-ARCHITECTURE.md
**Purpose:** Technical architecture design

**Contains:**
- System architecture diagram
- Technology stack decisions
- Database architecture
- API design
- Security architecture
- Scalability approach

---

#### 45-DATA-MODEL.md
**Purpose:** Core data model specification

**Contains:**
- Entity relationship diagram
- Entity definitions
- Relationship definitions
- Data validation rules
- Migration strategy

---

### Business (50-series)

#### 50-BUSINESS-MODEL.md
**Purpose:** Business model summary

**Contains:**
- Business model canvas
- Value creation logic
- Revenue model overview
- Key partnerships
- Cost structure overview

---

#### 51-REVENUE-STREAMS.md
**Purpose:** Revenue model detail

**Contains:**
- Revenue stream definitions
- Revenue mix targets
- Pricing rationale
- Revenue projections by stream

---

#### 52-PRICING.md
**Purpose:** Pricing strategy

**Contains:**
- Pricing philosophy
- Tier definitions
- Competitive pricing analysis
- Price sensitivity assessment
- Pricing evolution plan

---

#### 53-MARKETPLACE.md
**Purpose:** Marketplace dynamics

**Contains:**
- Marketplace model details
- Network effects analysis
- Chicken-and-egg strategy
- Quality control mechanisms
- Marketplace metrics

---

#### 54-UNIT-ECONOMICS.md
**Purpose:** Unit economics analysis

**Contains:**
- CAC by channel
- LTV by segment
- LTV:CAC ratios
- Payback period
- Contribution margin

---

### Go-to-Market (60-series)

#### 60-GTM-OVERVIEW.md
**Purpose:** Go-to-market summary

**Contains:**
- GTM strategy overview
- Target segment prioritization
- Channel strategy summary
- Launch approach
- Success metrics

---

#### 61-LAUNCH-STRATEGY.md
**Purpose:** Launch plan

**Contains:**
- Pre-launch activities
- Launch sequence
- Launch markets/segments
- Launch metrics
- Contingency plans

---

#### 62-CHANNELS.md
**Purpose:** Channel strategy

**Contains:**
- Channel inventory
- Channel prioritization
- Channel-specific tactics
- Channel economics
- Channel scaling plan

---

#### 63-MARKETING.md
**Purpose:** Marketing approach

**Contains:**
- Marketing strategy
- Messaging framework
- Content strategy
- Community building
- Marketing budget allocation

---

#### 64-PARTNERSHIPS.md
**Purpose:** Partnership strategy

**Contains:**
- Partnership categories
- Target partners
- Partnership value proposition
- Partnership economics
- Partnership roadmap

---

### Operations (70-series)

#### 70-OPERATIONS-OVERVIEW.md
**Purpose:** Operations summary

**Contains:**
- Operating model
- Key processes
- Technology operations
- Support model
- Compliance overview

---

#### 71-TEAM.md
**Purpose:** Team structure

**Contains:**
- Current team
- Hiring plan
- Organization structure
- Key roles and responsibilities
- Advisory board

---

#### 72-PROCESSES.md
**Purpose:** Key processes

**Contains:**
- Product development process
- Content curation process
- Customer support process
- Practitioner onboarding process
- Quality assurance process

---

#### 73-TECHNOLOGY.md
**Purpose:** Technology operations

**Contains:**
- Infrastructure
- Development practices
- Deployment process
- Monitoring and alerting
- Security operations

---

#### 74-COMPLIANCE.md
**Purpose:** Regulatory compliance

**Contains:**
- FDA software as medical device
- FTC health claims
- HIPAA requirements
- GDPR requirements
- Medical disclaimer approach

---

### Financials (80-series)

#### 80-FINANCIALS-OVERVIEW.md
**Purpose:** Financial summary

**Contains:**
- Financial highlights
- Key assumptions
- Funding status
- Use of funds
- Path to profitability

---

#### 81-PROJECTIONS.md
**Purpose:** Revenue projections

**Contains:**
- Revenue model assumptions
- Monthly/quarterly projections
- Scenario analysis (base/bull/bear)
- Sensitivity analysis

---

#### 82-COSTS.md
**Purpose:** Cost structure

**Contains:**
- Fixed costs
- Variable costs
- Cost drivers
- Cost optimization opportunities

---

#### 83-FUNDING.md
**Purpose:** Funding requirements

**Contains:**
- Current funding status
- Funding requirements by phase
- Use of funds breakdown
- Funding sources/strategy

---

#### 84-METRICS.md
**Purpose:** KPIs and metrics

**Contains:**
- North star metric
- User metrics
- Business metrics
- Product metrics
- Health outcome metrics

---

### Risk (90-series)

#### 90-RISK-OVERVIEW.md
**Purpose:** Risk summary

**Contains:**
- Risk framework
- Top 10 risks
- Risk mitigation summary
- Risk monitoring approach

---

#### 91-ASSUMPTIONS.md
**Purpose:** Critical assumptions

**Contains:**
- Market assumptions
- Product assumptions
- Business model assumptions
- Operational assumptions
- Assumption validation plan

---

#### 92-THREATS.md
**Purpose:** Threat analysis

**Contains:**
- Competitive threats
- Market threats
- Technology threats
- Regulatory threats
- Operational threats

---

#### 93-MITIGATIONS.md
**Purpose:** Risk mitigations

**Contains:**
- Mitigation strategies by risk
- Contingency plans
- Early warning indicators
- Response playbooks

---

## Document Lifecycle

### Status Definitions

| Status | Definition |
|--------|------------|
| Placeholder | Structure only, no content |
| Draft | Initial content, under development |
| Review | Content complete, under review |
| Final | Approved, ready for use |
| Deprecated | Outdated, being replaced |

### Update Cadence

| Document Type | Review Frequency |
|---------------|------------------|
| Strategy (10-series) | Quarterly |
| Market (20-series) | Quarterly |
| Customers (30-series) | Quarterly |
| Product (40-series) | Monthly |
| Business (50-series) | Quarterly |
| GTM (60-series) | Monthly |
| Operations (70-series) | As needed |
| Financials (80-series) | Monthly |
| Risk (90-series) | Quarterly |

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | — | Initial framework |
