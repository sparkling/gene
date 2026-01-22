# Jobs-to-be-Done Framework

**Document ID:** 33-JOBS-TO-BE-DONE
**Status:** Final
**Owner:** Product
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

JTBD framework identifies 8 core functional jobs, 6 emotional jobs, and 4 social jobs that drive platform usage. Primary functional job: "Understand what my genetic variants mean." Primary emotional job: "Feel in control of my health." Critical insight: Users aren't buying genetics analysis—they're hiring us to end their search for answers. The platform must satisfy both functional needs (information) and emotional needs (validation, belonging).

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| Primary functional job | "Understand what my genetic data means" | Universal entry point | Jan 2026 |
| Primary emotional job | "Feel in control of my health" | Deepest motivation | Jan 2026 |
| Job prioritization | Functional → Emotional → Social | Build trust, then community | Jan 2026 |

---

## JTBD Framework Overview

### What is JTBD?

Jobs-to-be-Done theory focuses on understanding what "job" customers are trying to accomplish, rather than focusing on product features or demographics.

**Key Principle:** People don't buy products; they "hire" products to do a job.

### JTBD Types

| Job Type | Description | Question |
|----------|-------------|----------|
| **Functional** | Practical tasks to accomplish | "What are you trying to get done?" |
| **Emotional** | How they want to feel | "How do you want to feel?" |
| **Social** | How they want to be perceived | "How do you want others to see you?" |

---

## Core Functional Jobs

### Job Map

```
                    FUNCTIONAL JOBS
                          │
    ┌─────────────────────┼─────────────────────┐
    │                     │                     │
UNDERSTAND            FIND                  TAKE
    │                     │                   ACTION
    │                     │                     │
┌───┴───┐             ┌───┴───┐            ┌───┴───┐
│Genetic│             │ Root  │            │Get rec│
│meaning│             │causes │            │ommen- │
│       │             │       │            │dations│
├───────┤             ├───────┤            ├───────┤
│Pathway│             │Connect│            │Track  │
│connec-│             │with   │            │progress│
│tions  │             │practi-│            │       │
│       │             │tioners│            │       │
└───────┘             └───────┘            └───────┘
```

### Functional Jobs Detail

| # | Job Statement | Priority | Feature Mapping |
|---|---------------|----------|-----------------|
| F1 | Understand what my genetic variants mean | P0 | SNP analysis, gene reports |
| F2 | Find root causes of my health issues | P0 | Pathway visualization, AI Q&A |
| F3 | Get actionable recommendations I can try | P0 | TCM/Ayurveda/nutrition suggestions |
| F4 | Learn about biochemistry relevant to me | P1 | Interactive pathways, education |
| F5 | Find practitioners who understand my situation | P1 | Practitioner marketplace |
| F6 | Connect with others who have similar issues | P1 | Community features |
| F7 | Track what works and what doesn't | P2 | Progress tracking, treatment history |
| F8 | Share findings with my healthcare team | P2 | PDF export, practitioner sharing |

### Functional Job Deep Dives

#### F1: Understand What My Genetic Variants Mean

| Dimension | Detail |
|-----------|--------|
| **Trigger** | Received 23andMe/AncestryDNA results, overwhelmed by raw data |
| **Current Solution** | Promethease ($12, overwhelming), Googling individual SNPs |
| **Pain Points** | Information overload, no guidance, conflicting sources |
| **Success Criteria** | Clear explanation of significant variants, prioritized insights |
| **Platform Solution** | AI-powered reports, plain language explanations, evidence links |

#### F2: Find Root Causes of My Health Issues

| Dimension | Detail |
|-----------|--------|
| **Trigger** | Years of symptoms, no diagnosis, dismissed by doctors |
| **Current Solution** | Self-research, multiple specialists, forums |
| **Pain Points** | Fragmented information, expensive, time-consuming |
| **Success Criteria** | Identify pathways/genes connected to symptoms |
| **Platform Solution** | Pathway visualization, symptom correlation, AI analysis |

#### F3: Get Actionable Recommendations I Can Try

| Dimension | Detail |
|-----------|--------|
| **Trigger** | Understand the problem, need solutions |
| **Current Solution** | Generic supplement advice, practitioner visits ($200+) |
| **Pain Points** | One-size-fits-all, expensive, contradictory |
| **Success Criteria** | Personalized, evidence-based, affordable |
| **Platform Solution** | TCM/Ayurveda herbs, nutritional cofactors, evidence-linked |

---

## Core Emotional Jobs

### Emotional Job Map

```
                    EMOTIONAL JOBS
                          │
    ┌─────────────────────┼─────────────────────┐
    │                     │                     │
VALIDATION            CONTROL               HOPE
    │                     │                     │
┌───┴───┐             ┌───┴───┐            ┌───┴───┐
│Symptoms│            │Agency │            │Answers│
│are real│            │over   │            │exist  │
│       │             │health │            │       │
├───────┤             ├───────┤            ├───────┤
│Not    │             │Under- │            │Progress│
│crazy  │             │stand  │            │is pos-│
│       │             │body   │            │sible  │
└───────┘             └───────┘            └───────┘
```

### Emotional Jobs Detail

| # | Job Statement | Intensity | Platform Response |
|---|---------------|-----------|-------------------|
| E1 | Feel in control of my health | HIGH | Self-service research, track progress |
| E2 | Feel validated that my symptoms are real | HIGH | Genetic evidence, pathway connections |
| E3 | Feel hopeful that answers exist | HIGH | Success stories, treatment options |
| E4 | Feel empowered with knowledge | MEDIUM | Educational content, visualizations |
| E5 | Feel less alone in my struggle | MEDIUM | Community, peer stories |
| E6 | Feel confident discussing with doctors | MEDIUM | PDF reports, evidence citations |

### Emotional Job Deep Dives

#### E1: Feel in Control of My Health

| Dimension | Detail |
|-----------|--------|
| **Underlying Need** | Agency, autonomy, empowerment |
| **Current State** | Dependent on medical system that has failed them |
| **Desired State** | Knowledgeable, capable of making informed decisions |
| **Platform Response** | Self-service tools, comprehensive knowledge access |
| **Risk** | Information without guidance can overwhelm |

#### E2: Feel Validated That My Symptoms Are Real

| Dimension | Detail |
|-----------|--------|
| **Underlying Need** | Recognition, credibility, self-worth |
| **Current State** | Dismissed by doctors ("It's all in your head") |
| **Desired State** | Genetic/biochemical evidence explaining symptoms |
| **Platform Response** | Gene-pathway-symptom connections, evidence citations |
| **Risk** | Over-attribution of symptoms to genetics |

---

## Core Social Jobs

### Social Job Map

```
                    SOCIAL JOBS
                          │
    ┌─────────────────────┼─────────────────────┐
    │                     │                     │
BELONGING             CREDIBILITY          CONTRIBUTION
    │                     │                     │
┌───┴───┐             ┌───┴───┐            ┌───┴───┐
│Connect│             │Inform-│            │Help   │
│with   │             │ed pat-│            │others │
│peers  │             │ient   │            │       │
├───────┤             ├───────┤            ├───────┤
│Part of│             │Health │            │Share  │
│commun-│             │advoca-│            │discov-│
│ity    │             │te     │            │eries  │
└───────┘             └───────┘            └───────┘
```

### Social Jobs Detail

| # | Job Statement | Priority | Platform Response |
|---|---------------|----------|-------------------|
| S1 | Connect with others who understand | P1 | Condition-specific communities |
| S2 | Be seen as an informed patient | P2 | Shareable reports, evidence backing |
| S3 | Help others who are struggling | P2 | Protocol sharing, community contributions |
| S4 | Contribute to collective knowledge | P3 | Anonymous data sharing, community protocols |

---

## Job-to-Feature Mapping

### Priority Matrix

| Job | Feature | Priority | Phase |
|-----|---------|----------|-------|
| F1: Understand genetics | SNP analysis, gene reports | P0 | MVP |
| F2: Find root causes | Pathway visualization | P0 | MVP |
| F3: Get recommendations | TCM/Ayurveda integration | P0 | MVP |
| E1: Feel in control | Self-service research | P0 | MVP |
| E2: Feel validated | Evidence-linked insights | P0 | MVP |
| F4: Learn biochemistry | Interactive pathways | P1 | Phase 2 |
| F5: Find practitioners | Practitioner marketplace | P1 | Phase 2 |
| S1: Connect with peers | Community features | P1 | Phase 2 |
| E5: Feel less alone | Peer communities | P1 | Phase 2 |
| F6: Connect with others | Geographic communities | P2 | Phase 3 |
| F7: Track progress | Treatment tracking | P2 | Phase 3 |
| F8: Share with doctors | Export, sharing | P2 | Phase 3 |

### Feature Alignment Check

| Feature | Jobs Served | Alignment Score |
|---------|-------------|-----------------|
| AI-powered Q&A | F1, F2, F3, F4 | HIGH |
| Pathway visualization | F1, F2, F4 | HIGH |
| TCM herb database | F3 | HIGH |
| Ayurveda integration | F3 | HIGH |
| Community features | F6, S1, S3, E5 | HIGH |
| Practitioner marketplace | F5 | MEDIUM |
| Progress tracking | F7 | MEDIUM |
| PDF export | F8 | LOW |

---

## User Journey by Job

### The "Unhelped" Journey

```
STAGE 1: AWARENESS              STAGE 2: EXPLORATION
"Something is wrong"            "I need to understand"
        │                               │
   ┌────┴────┐                    ┌─────┴─────┐
   │ Symptoms│                    │  Upload   │
   │ persist │                    │   DNA     │
   │         │                    │   data    │
   └────┬────┘                    └─────┬─────┘
        │                               │
        ▼                               ▼
JOBS: E2 (validation)           JOBS: F1 (understand)
      E3 (hope)                       F2 (root causes)

STAGE 3: LEARNING               STAGE 4: ACTION
"Now I see connections"         "I can try something"
        │                               │
   ┌────┴────┐                    ┌─────┴─────┐
   │ Explore │                    │ Get recs  │
   │pathways │                    │ Find help │
   │         │                    │           │
   └────┬────┘                    └─────┬─────┘
        │                               │
        ▼                               ▼
JOBS: F4 (learn)                JOBS: F3 (recommendations)
      E1 (control)                    F5 (practitioners)
      E4 (empowered)                  S1 (community)

STAGE 5: COMMUNITY
"I'm not alone"
        │
   ┌────┴────┐
   │ Connect │
   │  Share  │
   │ Contrib │
   └────┬────┘
        │
        ▼
JOBS: F6 (connect)
      S3 (help others)
      E5 (belonging)
```

---

## Outcome Expectations

### Functional Outcomes

| Job | Expected Outcome | Measurement |
|-----|------------------|-------------|
| F1 | Understand 80%+ of significant variants | Report comprehension survey |
| F2 | Identify at least 1 plausible root cause | User-reported insights |
| F3 | Receive 5+ actionable recommendations | Recommendation count |
| F4 | Learn pathway mechanics | Education completion |
| F5 | Find 3+ relevant practitioners | Search success rate |

### Emotional Outcomes

| Job | Expected Outcome | Measurement |
|-----|------------------|-------------|
| E1 | Feel informed and capable | Post-use survey |
| E2 | Feel symptoms are validated | Sentiment analysis |
| E3 | Feel hopeful about progress | NPS score |
| E5 | Feel connected to community | Engagement metrics |

### Social Outcomes

| Job | Expected Outcome | Measurement |
|-----|------------------|-------------|
| S1 | Find peer connections | Community joins |
| S3 | Contribute to community | Protocol shares |

---

## Competitive JTBD Analysis

### How Competitors Address Jobs

| Job | SelfDecode | StrateGene | LiveWello | Gene Platform |
|-----|------------|------------|-----------|---------------|
| F1: Understand genetics | ✅ | ✅ | ⚠️ | ✅ |
| F2: Find root causes | ⚠️ | ✅ | ❌ | ✅ |
| F3: Get recommendations | ⚠️ Supplements | ❌ | ❌ | ✅ TCM/Ayur |
| F5: Find practitioners | ❌ | ❌ | ❌ | ✅ |
| E2: Feel validated | ⚠️ | ✅ | ❌ | ✅ |
| S1: Connect with peers | ❌ | ❌ | ⚠️ | ✅ |

**Gap:** No competitor fully addresses practitioner finding (F5) or community connection (S1).

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [31-SEGMENTS](./31-SEGMENTS.md) | Segment definitions |
| [32-PERSONAS](./32-PERSONAS.md) | Persona detail |
| [34-USER-JOURNEYS](./34-USER-JOURNEYS.md) | Journey mapping |
| [41-FEATURES](../40-product/41-features.md) | Feature alignment |

---

## Open Questions

- [ ] Validate job prioritization with user research
- [ ] Quantify emotional job intensity across segments
- [ ] Map social jobs to community features

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Product | Complete JTBD framework |
