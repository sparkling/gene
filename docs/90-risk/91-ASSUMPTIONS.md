# Critical Assumptions

**Document ID:** 91-ASSUMPTIONS
**Status:** Final
**Owner:** Strategy
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Platform success depends on 24 critical assumptions across 5 categories: Market, Product, Business Model, Technical, and Operational. Highest-risk assumptions: user willingness to pay for genetics + traditional medicine integration (CRITICAL), TCM/Ayurveda data quality sufficient for recommendations (HIGH), practitioner marketplace achieves liquidity (HIGH). Validation plan prioritizes these through beta testing and early user research.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Assumption framework | 5-category model | MECE coverage of business drivers | Jan 2026 |
| Validation priority | Market assumptions first | Must validate demand before building | Jan 2026 |
| Confidence threshold | 70%+ to proceed | Balance speed with risk | Jan 2026 |
| Review cadence | Monthly | Fast iteration on learnings | Jan 2026 |

---

## Assumption Framework

```
                          ASSUMPTION CATEGORIES
                                  │
    ┌────────────┬────────────────┼────────────────┬────────────┐
    │            │                │                │            │
  MARKET      PRODUCT        BUSINESS         TECHNICAL    OPERATIONAL
    │            │            MODEL               │            │
    │            │                │                │            │
 Demand       Features         Pricing          Data       Team
 Segments     Adoption         Revenue          Quality    Scaling
 Timing       Value            Take Rate        AI         Support
 Competition  Differentiation  Conversion       Integration
```

---

## Assumption Registry

### Category 1: Market Assumptions

| ID | Assumption | Confidence | Risk if Wrong | Validation Method |
|----|------------|------------|---------------|-------------------|
| M1 | Users will pay for genetics + traditional medicine | 60% | CRITICAL | Beta pricing test |
| M2 | Chronic illness segment is large and underserved | 80% | HIGH | User research, forum analysis |
| M3 | Traditional medicine interest is growing | 75% | MEDIUM | Trend data, search volume |
| M4 | Competitors won't add TCM/Ayurveda in <12mo | 50% | HIGH | Competitive monitoring |
| M5 | DTC genetic testing market continues growing | 85% | LOW | Industry reports |
| M6 | Users trust AI for health recommendations | 65% | MEDIUM | User interviews |

### Category 2: Product Assumptions

| ID | Assumption | Confidence | Risk if Wrong | Validation Method |
|----|------------|------------|---------------|-------------------|
| P1 | Interactive pathways are valued over static PDFs | 70% | MEDIUM | A/B testing |
| P2 | AI Q&A delivers accurate, useful responses | 65% | HIGH | Accuracy testing, user feedback |
| P3 | Gene-herb-condition mapping is core differentiator | 80% | HIGH | Competitive analysis |
| P4 | Users want community features | 60% | MEDIUM | Feature request tracking |
| P5 | 200 SNPs sufficient for MVP value | 75% | MEDIUM | Expert review, user feedback |
| P6 | Cytoscape.js suitable for pathway visualization | 80% | LOW | Technical prototype |

### Category 3: Business Model Assumptions

| ID | Assumption | Confidence | Risk if Wrong | Validation Method |
|----|------------|------------|---------------|-------------------|
| B1 | Freemium converts at 10%+ | 55% | HIGH | Beta conversion tracking |
| B2 | $79/year is acceptable price point | 60% | HIGH | Price sensitivity testing |
| B3 | Practitioner marketplace achieves liquidity | 50% | HIGH | Early practitioner outreach |
| B4 | 15-20% take rate acceptable to practitioners | 70% | MEDIUM | Practitioner interviews |
| B5 | Report template marketplace has demand | 55% | MEDIUM | Creator interest survey |
| B6 | LTV:CAC achieves 3:1+ | 60% | HIGH | Unit economics tracking |

### Category 4: Technical Assumptions

| ID | Assumption | Confidence | Risk if Wrong | Validation Method |
|----|------------|------------|---------------|-------------------|
| T1 | TCM database quality sufficient for recommendations | 55% | HIGH | Data validation, expert review |
| T2 | DNA parsing handles major formats reliably | 75% | MEDIUM | Test with real files |
| T3 | AI query costs manageable at scale | 60% | HIGH | Cost modeling, caching tests |
| T4 | Neo4j handles relationship queries efficiently | 80% | LOW | Load testing |
| T5 | RAG architecture delivers quality answers | 65% | HIGH | Response quality testing |
| T6 | Vector search scales to 1M+ users | 75% | LOW | Benchmark testing |

### Category 5: Operational Assumptions

| ID | Assumption | Confidence | Risk if Wrong | Validation Method |
|----|------------|------------|---------------|-------------------|
| O1 | Can hire 2 full-stack engineers in 30 days | 60% | MEDIUM | Recruiting pipeline |
| O2 | Solo founder can manage Phase 1 operations | 65% | MEDIUM | Time tracking, delegation |
| O3 | Customer support scalable via community | 55% | MEDIUM | Community engagement metrics |
| O4 | Medical/content expert available part-time | 70% | MEDIUM | Network outreach |
| O5 | Legal/compliance reviewable via fractional counsel | 75% | LOW | Legal firm engagement |
| O6 | Can launch MVP in 4 months | 55% | HIGH | Milestone tracking |

---

## Critical Assumptions (Must Validate First)

### Assumption M1: Users Will Pay for Genetics + Traditional Medicine

**Statement:** Users will pay $79-179/year for a platform connecting their genetic data to traditional medicine recommendations.

| Factor | Assessment |
|--------|------------|
| **Confidence** | 60% |
| **Risk if Wrong** | CRITICAL — No business |
| **Supporting Evidence** | SelfDecode ($199-899), StrateGene ($95), ADNTRO ($200-400) prove willingness to pay for genetic analysis |
| **Concerning Evidence** | TCM/Ayurveda angle unvalidated; may be niche |
| **Validation Plan** | Beta pricing test with 500+ users |
| **Decision Point** | If <5% conversion at $79, reconsider pricing or positioning |

### Assumption T1: TCM Database Quality Sufficient

**Statement:** Available TCM databases (TCMSP, HERB, BATMAN-TCM) have sufficient quality for accurate gene-herb recommendations.

| Factor | Assessment |
|--------|------------|
| **Confidence** | 55% |
| **Risk if Wrong** | HIGH — Core feature compromised |
| **Supporting Evidence** | Academic databases exist with thousands of herb-target connections |
| **Concerning Evidence** | Data may be incomplete, outdated, or inconsistent |
| **Validation Plan** | Expert review of top 100 herbs; cross-reference multiple sources |
| **Decision Point** | If >20% of top herbs have unreliable data, narrow initial scope |

### Assumption B3: Practitioner Marketplace Achieves Liquidity

**Statement:** The practitioner marketplace will achieve sufficient supply (50+ practitioners) and demand (200+ bookings/month) within 6 months of launch.

| Factor | Assessment |
|--------|------------|
| **Confidence** | 50% |
| **Risk if Wrong** | HIGH — Marketplace revenue stream fails |
| **Supporting Evidence** | Practitioners need client acquisition; platform offers unique value |
| **Concerning Evidence** | Chicken-egg problem; practitioners won't join without users |
| **Validation Plan** | Pre-seed 20 practitioners before user launch; offer incentives |
| **Decision Point** | If <10 practitioners committed by launch, delay marketplace |

---

## Assumption Confidence Matrix

```
                    CONFIDENCE
                    Low    Med    High
                    (<60%) (60-75%) (>75%)
              ┌─────────┬─────────┬─────────┐
   CRITICAL   │ M1, B3  │         │         │
              ├─────────┼─────────┼─────────┤
IMPACT HIGH   │ T1, T3  │ P2, B1  │ M2, P3  │
              │ B5, O6  │ B2, B6  │         │
              ├─────────┼─────────┼─────────┤
   MEDIUM     │ P4, O3  │ M4, M6  │ P1, P5  │
              │         │ T5, O1  │ B4, T2  │
              │         │ O2, O4  │         │
              ├─────────┼─────────┼─────────┤
   LOW        │         │         │ M3, M5  │
              │         │         │ P6, T4  │
              │         │         │ T6, O5  │
              └─────────┴─────────┴─────────┘

PRIORITY: Validate top-left quadrant FIRST
```

---

## Validation Priority Stack

### Tier 1: Validate Before Building (Pre-MVP)

| ID | Assumption | Validation Method | Timeline | Owner |
|----|------------|-------------------|----------|-------|
| M1 | Willingness to pay | Landing page test, user interviews | Week 1-2 | Founder |
| M2 | Segment underserved | Forum analysis, user interviews | Week 1-2 | Founder |
| T1 | TCM data quality | Expert review of top 50 herbs | Week 2-4 | Data |

### Tier 2: Validate During MVP (Phase 1)

| ID | Assumption | Validation Method | Timeline | Owner |
|----|------------|-------------------|----------|-------|
| P2 | AI accuracy | User feedback, accuracy metrics | Month 2-4 | Engineering |
| T3 | AI costs | Cost tracking, caching impact | Month 2-4 | Engineering |
| B1 | Freemium conversion | Beta conversion tracking | Month 4+ | Product |
| B2 | Price point | Price testing | Month 4+ | Product |
| O6 | 4-month MVP | Milestone tracking | Month 1-4 | Founder |

### Tier 3: Validate Post-MVP (Phase 2+)

| ID | Assumption | Validation Method | Timeline | Owner |
|----|------------|-------------------|----------|-------|
| B3 | Marketplace liquidity | Practitioner/booking metrics | Month 5-8 | Product |
| B6 | Unit economics | LTV:CAC calculation | Month 6+ | Finance |
| P4 | Community demand | Feature usage, requests | Month 6+ | Product |

---

## Assumption Testing Methods

### A/B Testing

| Assumption | Test Design | Success Metric |
|------------|-------------|----------------|
| P1 | Interactive vs static pathway | Time on page, satisfaction |
| B2 | $59 vs $79 vs $99 pricing | Conversion rate |
| B1 | Free tier feature limits | Upgrade rate |

### User Research

| Assumption | Research Method | Sample Size |
|------------|-----------------|-------------|
| M1 | Willingness-to-pay survey | 200+ |
| M2 | Problem interviews | 30+ |
| M6 | Trust in AI study | 50+ |

### Expert Review

| Assumption | Expert Type | Review Scope |
|------------|-------------|--------------|
| T1 | TCM practitioner | Top 100 herbs |
| T1 | Data scientist | Data completeness |
| P5 | Geneticist | SNP panel adequacy |

### Quantitative Tracking

| Assumption | Metric | Target |
|------------|--------|--------|
| B1 | Conversion rate | >10% |
| B6 | LTV:CAC | >3:1 |
| T3 | Cost per AI query | <$0.05 |

---

## Assumption Update Log

| Date | Assumption | Old Confidence | New Confidence | Reason |
|------|------------|----------------|----------------|--------|
| — | — | — | — | Initial baseline |

*Updated monthly based on validation results*

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [90-RISK-OVERVIEW](./90-RISK-OVERVIEW.md) | Risks depend on assumptions |
| [92-THREATS](./92-THREATS.md) | Threat assumptions |
| [42-ROADMAP](../40-product/42-ROADMAP.md) | Timeline assumptions |
| [54-UNIT-ECONOMICS](../50-business/54-UNIT-ECONOMICS.md) | Financial assumptions |

---

## Open Questions

- [ ] How to validate TCM efficacy claims without clinical trials?
- [ ] What minimum practitioner count for marketplace viability?
- [ ] Acceptable error rate for AI health recommendations?
- [ ] International market assumptions for Phase 3+?

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Strategy | Complete assumptions registry |
