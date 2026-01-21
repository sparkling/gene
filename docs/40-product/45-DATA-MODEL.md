# Core Data Model

**Document ID:** 45-DATA-MODEL
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

**Three-database architecture** using PostgreSQL for HIPAA-compliant user data, RuVector Graph for knowledge graph relationships, and RuVector PostgreSQL for vector embeddings. Core entities: User, GeneticProfile, SNP, Gene, Pathway, Nutrient, Herb, Condition, Practitioner. The knowledge graph connects THREE WORLDS: genes/SNPs (genetics) - herbs/formulas (traditional medicine) - nutrients/foods (nutrition). RuVector Graph provides Cypher query compatibility and hypergraph support for complex multi-node relationships. Designed for HIPAA compliance with field-level encryption on sensitive genetic data.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| User data storage | PostgreSQL (port 5433) | ACID, relational integrity, HIPAA isolation | Jan 2026 |
| Knowledge relationships | RuVector Graph (@ruvector/graph-node) | Unified stack, Cypher compatible, 131K ops/sec, hypergraphs | Jan 2026 |
| Vector embeddings | RuVector PostgreSQL (port 5432) | HNSW indexing, 150x-12,500x faster search | Jan 2026 |
| Genetic data encryption | Field-level AES-256 | HIPAA requirement | Jan 2026 |
| ID strategy | UUIDs | Distribution-friendly, privacy | Jan 2026 |

---

## Data Model Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                 USER DOMAIN (PostgreSQL - Port 5433)             │
│                      HIPAA Compliant, Encrypted                  │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  User ──┬── GeneticProfile[] ──── UserSNP[]                     │
│         ├── LabResult[]                                          │
│         ├── Symptom[]                                            │
│         ├── Treatment[]                                          │
│         ├── Subscription                                         │
│         └── Order[]                                              │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ Links via rsid (SNP identifier)
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│              KNOWLEDGE DOMAIN (RuVector Graph)                   │
│              @ruvector/graph-node - 131K ops/sec                 │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Gene ──── SNP ──── Pathway ──── Nutrient                       │
│    │         │         │            │                            │
│    └─────────┴─────────┴────────────┘                           │
│                  │                                               │
│              Condition ──── Herb ──── Formula                    │
│                  │            │          │                       │
│              Publication  Modality   Modality                    │
│                                                                  │
│  ═══════════════════════════════════════════════════════════    │
│  HYPEREDGES (n-ary relationships):                               │
│  • METHYLATION_CYCLE_PARTICIPANTS: [MTHFR, MTR, MTRR, SHMT1]    │
│  • FORMULA_COMPOSITION: [huang_qi, dang_shen, bai_zhu, gan_cao] │
│  ═══════════════════════════════════════════════════════════    │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
                              │
                              │ Vector embeddings for RAG
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│              VECTOR DOMAIN (RuVector PostgreSQL - Port 5432)     │
│              HNSW Indexing - 150x-12,500x faster search          │
├─────────────────────────────────────────────────────────────────┤
│  • embeddings (research papers, clinical notes)                  │
│  • patterns (learned AI patterns)                                │
│  • hyperbolic_embeddings (hierarchical data)                     │
└─────────────────────────────────────────────────────────────────┘
```

---

## PostgreSQL Schema (User Domain)

### User Entity

```sql
CREATE TABLE users (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    email VARCHAR(255) UNIQUE NOT NULL,
    email_verified BOOLEAN DEFAULT FALSE,
    password_hash VARCHAR(255),

    -- Profile
    first_name VARCHAR(100),
    last_name VARCHAR(100),
    date_of_birth DATE,
    biological_sex VARCHAR(20), -- 'male', 'female', 'other'
    timezone VARCHAR(50) DEFAULT 'UTC',

    -- Account
    role VARCHAR(20) DEFAULT 'user', -- 'user', 'practitioner', 'admin'
    subscription_tier VARCHAR(20) DEFAULT 'free',
    subscription_status VARCHAR(20) DEFAULT 'active',

    -- Preferences
    preferences JSONB DEFAULT '{}',

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    last_login_at TIMESTAMPTZ,
    deleted_at TIMESTAMPTZ -- soft delete
);

CREATE INDEX idx_users_email ON users(email);
CREATE INDEX idx_users_role ON users(role);
```

### GeneticProfile Entity

```sql
CREATE TABLE genetic_profiles (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Source
    source VARCHAR(50) NOT NULL, -- '23andme', 'ancestry', 'vcf', 'manual'
    source_version VARCHAR(20),
    original_filename VARCHAR(255),

    -- Processing
    status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'processing', 'complete', 'failed'
    snp_count INTEGER,
    processing_started_at TIMESTAMPTZ,
    processing_completed_at TIMESTAMPTZ,
    error_message TEXT,

    -- Encrypted raw data (HIPAA)
    raw_data_encrypted BYTEA, -- AES-256 encrypted
    raw_data_iv BYTEA,

    -- Metadata
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),

    CONSTRAINT fk_user FOREIGN KEY (user_id) REFERENCES users(id)
);

CREATE INDEX idx_genetic_profiles_user ON genetic_profiles(user_id);
CREATE INDEX idx_genetic_profiles_status ON genetic_profiles(status);
```

### UserSNP Entity

```sql
CREATE TABLE user_snps (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    genetic_profile_id UUID NOT NULL,

    -- SNP Data
    rsid VARCHAR(20) NOT NULL, -- e.g., 'rs1801133'
    chromosome VARCHAR(5),
    position BIGINT,
    genotype VARCHAR(10), -- e.g., 'CT', 'AA'

    -- Reference (links to Neo4j)
    dbsnp_id VARCHAR(20),

    -- Interpretation
    risk_level VARCHAR(20), -- 'normal', 'heterozygous', 'homozygous_risk'
    clinical_significance VARCHAR(50),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),

    CONSTRAINT fk_genetic_profile FOREIGN KEY (genetic_profile_id)
        REFERENCES genetic_profiles(id) ON DELETE CASCADE,
    CONSTRAINT unique_profile_rsid UNIQUE (genetic_profile_id, rsid)
);

CREATE INDEX idx_user_snps_profile ON user_snps(genetic_profile_id);
CREATE INDEX idx_user_snps_rsid ON user_snps(rsid);
CREATE INDEX idx_user_snps_risk ON user_snps(risk_level);
```

### LabResult Entity

```sql
CREATE TABLE lab_results (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Test Info
    test_name VARCHAR(255) NOT NULL,
    test_code VARCHAR(50),
    category VARCHAR(50), -- 'blood', 'hormones', 'stool', 'urine'

    -- Result
    value DECIMAL(10, 4),
    value_text VARCHAR(255), -- for non-numeric results
    unit VARCHAR(50),
    reference_range_low DECIMAL(10, 4),
    reference_range_high DECIMAL(10, 4),

    -- Interpretation
    status VARCHAR(20), -- 'normal', 'low', 'high', 'critical'

    -- Metadata
    test_date DATE,
    lab_name VARCHAR(255),
    notes TEXT,

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_lab_results_user ON lab_results(user_id);
CREATE INDEX idx_lab_results_category ON lab_results(category);
CREATE INDEX idx_lab_results_date ON lab_results(test_date);
```

### Symptom Entity

```sql
CREATE TABLE symptoms (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Symptom
    name VARCHAR(255) NOT NULL,
    category VARCHAR(50), -- 'physical', 'mental', 'digestive', 'neurological'

    -- Severity & Tracking
    severity INTEGER CHECK (severity BETWEEN 1 AND 10),
    frequency VARCHAR(50), -- 'constant', 'daily', 'weekly', 'occasional'

    -- Timing
    started_at DATE,
    ended_at DATE,

    -- Notes
    notes TEXT,
    triggers TEXT[],

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_symptoms_user ON symptoms(user_id);
CREATE INDEX idx_symptoms_category ON symptoms(category);
```

### Treatment Entity

```sql
CREATE TABLE treatments (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Treatment
    name VARCHAR(255) NOT NULL,
    type VARCHAR(50), -- 'supplement', 'medication', 'herb', 'lifestyle'
    modality VARCHAR(50), -- 'western', 'tcm', 'ayurveda', 'kampo'

    -- Dosage
    dosage VARCHAR(100),
    frequency VARCHAR(100),

    -- Duration
    started_at DATE,
    ended_at DATE,

    -- Effectiveness
    effectiveness_rating INTEGER CHECK (effectiveness_rating BETWEEN 1 AND 5),
    side_effects TEXT[],
    notes TEXT,

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_treatments_user ON treatments(user_id);
CREATE INDEX idx_treatments_type ON treatments(type);
CREATE INDEX idx_treatments_modality ON treatments(modality);
```

### Subscription Entity

```sql
CREATE TABLE subscriptions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Plan
    tier VARCHAR(20) NOT NULL, -- 'free', 'report', 'annual', 'pro'
    price_cents INTEGER,
    billing_period VARCHAR(20), -- 'one_time', 'monthly', 'yearly'

    -- Status
    status VARCHAR(20) DEFAULT 'active', -- 'active', 'canceled', 'past_due', 'expired'

    -- Stripe
    stripe_subscription_id VARCHAR(255),
    stripe_customer_id VARCHAR(255),

    -- Dates
    current_period_start TIMESTAMPTZ,
    current_period_end TIMESTAMPTZ,
    canceled_at TIMESTAMPTZ,

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),

    CONSTRAINT unique_user_subscription UNIQUE (user_id)
);

CREATE INDEX idx_subscriptions_user ON subscriptions(user_id);
CREATE INDEX idx_subscriptions_status ON subscriptions(status);
CREATE INDEX idx_subscriptions_stripe ON subscriptions(stripe_subscription_id);
```

### Order Entity (Marketplace)

```sql
CREATE TABLE orders (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id),

    -- Type
    order_type VARCHAR(50), -- 'subscription', 'report_template', 'consultation'

    -- Items
    items JSONB NOT NULL, -- array of {product_id, name, quantity, price_cents}

    -- Pricing
    subtotal_cents INTEGER NOT NULL,
    discount_cents INTEGER DEFAULT 0,
    tax_cents INTEGER DEFAULT 0,
    total_cents INTEGER NOT NULL,

    -- Payment
    status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'paid', 'failed', 'refunded'
    stripe_payment_intent_id VARCHAR(255),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    paid_at TIMESTAMPTZ,
    refunded_at TIMESTAMPTZ
);

CREATE INDEX idx_orders_user ON orders(user_id);
CREATE INDEX idx_orders_status ON orders(status);
CREATE INDEX idx_orders_created ON orders(created_at);
```

### Practitioner Entity

```sql
CREATE TABLE practitioners (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES users(id) ON DELETE CASCADE,

    -- Profile
    display_name VARCHAR(255) NOT NULL,
    bio TEXT,
    profile_image_url VARCHAR(500),

    -- Credentials
    credentials TEXT[], -- ['MD', 'ND', 'LAc']
    license_number VARCHAR(100),
    license_state VARCHAR(50),
    license_verified BOOLEAN DEFAULT FALSE,

    -- Specialties
    specialties TEXT[], -- ['functional_medicine', 'tcm', 'ayurveda']
    conditions TEXT[], -- ['me_cfs', 'autoimmune', 'adhd']
    modalities TEXT[], -- ['tcm', 'ayurveda', 'western_herbal']

    -- Availability
    accepting_new_clients BOOLEAN DEFAULT TRUE,
    consultation_rate_cents INTEGER,
    availability_notes TEXT,

    -- Marketplace
    total_consultations INTEGER DEFAULT 0,
    average_rating DECIMAL(3, 2),
    review_count INTEGER DEFAULT 0,

    -- Status
    status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'approved', 'suspended'

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),

    CONSTRAINT unique_practitioner_user UNIQUE (user_id)
);

CREATE INDEX idx_practitioners_status ON practitioners(status);
CREATE INDEX idx_practitioners_specialties ON practitioners USING GIN(specialties);
CREATE INDEX idx_practitioners_modalities ON practitioners USING GIN(modalities);
```

---

## RuVector Graph Schema (Knowledge Domain)

**Note:** RuVector Graph uses `@ruvector/graph-node` which provides Cypher-compatible syntax. The schema below uses the same query language as Neo4j for easy migration.

### Node Types

#### Gene Node

```typescript
// Using RuVector Graph API
import { CodeGraph } from 'ruvector/dist/core/graph-wrapper';

const graph = new CodeGraph({ storagePath: './data/knowledge-graph.db' });

graph.createNode('HGNC:7436', ['Gene'], {
    symbol: 'MTHFR',
    name: 'Methylenetetrahydrofolate reductase',
    chromosome: '1',
    description: 'Enzyme in folate metabolism',
    ensembl_id: 'ENSG00000177000',
    uniprot_id: 'P42898'
});
```

#### SNP Node

```typescript
graph.createNode('rs1801133', ['SNP'], {
    rsid: 'rs1801133',
    gene_symbol: 'MTHFR',
    chromosome: '1',
    position: 11856378,
    ref_allele: 'G',
    alt_allele: 'A',
    minor_allele: 'A',
    maf: 0.34,
    clinical_significance: 'drug_response',
    common_name: 'MTHFR C677T'
});
```

#### Pathway Node

```typescript
graph.createNode('REACT_R-HSA-196741', ['Pathway'], {
    name: 'Methylation Cycle',
    description: 'One-carbon metabolism and methyl group transfer',
    category: 'metabolism',
    reactome_id: 'R-HSA-196741',
    kegg_id: 'hsa00670'
});
```

#### Nutrient Node

```typescript
graph.createNode('NUT_methylfolate', ['Nutrient'], {
    name: 'Methylfolate',
    alternate_names: ['5-MTHF', 'L-methylfolate', 'Metafolin'],
    category: 'vitamin',
    subcategory: 'B_vitamin',
    usda_id: '1177',
    recommended_daily: '400-800mcg'
});
```

#### Herb Node

```typescript
graph.createNode('TCM_huangqi', ['Herb'], {
    name: 'Huang Qi',
    latin_name: 'Astragalus membranaceus',
    english_name: 'Astragalus',
    modality: 'tcm',
    tcmsp_id: 'MOL000239',
    properties: ['qi_tonic', 'immune_support'],
    meridians: ['spleen', 'lung'],
    temperature: 'warm',
    taste: 'sweet'
});
```

#### Condition Node

```typescript
graph.createNode('COND_me_cfs', ['Condition'], {
    name: 'ME/CFS',
    full_name: 'Myalgic Encephalomyelitis / Chronic Fatigue Syndrome',
    icd10_code: 'G93.32',
    mesh_id: 'D015673',
    category: 'neurological'
});
```

#### Publication Node

```typescript
graph.createNode('PMID_12345678', ['Publication'], {
    pmid: '12345678',
    title: 'MTHFR polymorphisms and folate status',
    authors: ['Smith J', 'Jones A'],
    journal: 'Nature Genetics',
    year: 2023,
    doi: '10.1038/ng.1234',
    evidence_level: 'clinical_trial'
});
```

#### Formula Node (TCM/Kampo)

```typescript
graph.createNode('FORMULA_buzhongyiqi', ['Formula'], {
    name: 'Bu Zhong Yi Qi Tang',
    english_name: 'Tonify the Middle and Augment the Qi Decoction',
    modality: 'tcm',
    category: 'qi_tonifying',
    indication: 'spleen_qi_deficiency'
});
```

### Relationship Types (Edges)

```typescript
// Gene → SNP
graph.createEdge('HGNC:7436', 'rs1801133', 'HAS_VARIANT');

// Gene → Pathway
graph.createEdge('HGNC:7436', 'REACT_R-HSA-196741', 'INVOLVED_IN', { role: 'enzyme' });

// SNP → Nutrient
graph.createEdge('rs1801133', 'NUT_methylfolate', 'AFFECTS', {
    effect: 'reduces_activity',
    magnitude: 'moderate'
});

// Pathway → Nutrient
graph.createEdge('REACT_R-HSA-196741', 'NUT_methylfolate', 'REQUIRES_COFACTOR');

// Herb → Condition
graph.createEdge('TCM_huangqi', 'COND_me_cfs', 'TREATS', {
    evidence_level: 'clinical',
    mechanism: 'anti_inflammatory'
});

// Herb → Formula
graph.createEdge('TCM_huangqi', 'FORMULA_buzhongyiqi', 'INGREDIENT_OF', { role: 'chief' });

// SNP → Condition
graph.createEdge('rs1801133', 'COND_depression', 'ASSOCIATED_WITH', {
    odds_ratio: 1.5,
    p_value: 0.001
});

// SNP/Herb → Publication
graph.createEdge('rs1801133', 'PMID_12345678', 'HAS_RESEARCH');
graph.createEdge('TCM_huangqi', 'PMID_12345678', 'HAS_RESEARCH');
```

### Hyperedges (Multi-Node Relationships)

RuVector Graph uniquely supports hyperedges - edges connecting multiple nodes simultaneously:

```typescript
// Multiple genes participating in a pathway
graph.createHyperedge(
    ['HGNC:7436', 'HGNC:7437', 'HGNC:7438', 'HGNC:6027'],  // MTHFR, MTR, MTRR, SHMT1
    'METHYLATION_CYCLE_PARTICIPANTS',
    {
        pathway_id: 'REACT_R-HSA-196741',
        role: 'core_enzymes',
        interaction_type: 'sequential'
    }
);

// TCM formula composition
graph.createHyperedge(
    ['TCM_huangqi', 'TCM_dangshen', 'TCM_baizhu', 'TCM_gancao'],
    'FORMULA_COMPOSITION',
    {
        formula: 'Bu Zhong Yi Qi Tang',
        traditional_use: 'qi_tonifying',
        chief_herb: 'TCM_huangqi'
    }
);
```

### Sample Queries (Cypher Compatible)

#### Find Nutrients Affected by User's SNPs

```typescript
const result = graph.cypher(`
    MATCH (snp:SNP)-[a:AFFECTS]->(nutrient:Nutrient)
    WHERE snp.rsid IN $user_rsids
    RETURN snp.rsid, snp.common_name, nutrient.name, a.effect, a.magnitude
    ORDER BY a.magnitude DESC
`, { user_rsids: ['rs1801133', 'rs1801131', 'rs1805087'] });
```

#### Find TCM Herbs for a Condition via Pathways

```typescript
const result = graph.cypher(`
    MATCH (condition:Condition {name: $condition_name})
          <-[:TREATS]-(herb:Herb {modality: 'tcm'})
          -[:HAS_RESEARCH]->(pub:Publication)
    WHERE pub.evidence_level IN ['clinical_trial', 'meta_analysis']
    RETURN herb.name, herb.properties, count(pub) as research_count
    ORDER BY research_count DESC
    LIMIT 10
`, { condition_name: 'ME/CFS' });
```

#### Find Pathway Genes with User's Risk Variants

```typescript
const result = graph.cypher(`
    MATCH (pathway:Pathway {name: $pathway_name})
          <-[:INVOLVED_IN]-(gene:Gene)
          -[:HAS_VARIANT]->(snp:SNP)
    WHERE snp.rsid IN $user_risk_snps
    RETURN pathway.name, gene.symbol, snp.rsid, snp.clinical_significance
`, { pathway_name: 'Methylation Cycle', user_risk_snps: ['rs1801133', 'rs1801131'] });
```

### Graph Algorithms

```typescript
// Find most connected genes (hub identification)
const pageRanks = graph.pageRank(20, 0.85);

// Detect communities of related entities
const communities = graph.communities();  // Louvain algorithm

// Find shortest path between gene and condition
const path = graph.shortestPath('HGNC:7436', 'COND_depression', 5);

// Get all mechanistic paths (for showing multiple pathways)
const allPaths = graph.allPaths('HGNC:7436', 'COND_depression', 4, 10);

// Find betweenness centrality (key connector nodes)
const centrality = graph.betweennessCentrality();
```

---

## Entity Relationship Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    COMPLETE DATA MODEL                                   │
│            Three-Database Architecture with RuVector                     │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  PostgreSQL (User Domain - Port 5433, HIPAA Compliant)                   │
│  ══════════════════════════════════════════════════                     │
│                                                                          │
│  User ─────────┬────────────────────────────────────────┐               │
│    │           │                                         │               │
│    │     GeneticProfile                            Subscription         │
│    │           │                                         │               │
│    │       UserSNP[]──────────────────────┬─────────────┼─┐             │
│    │           │                          │              │ │             │
│    │     LabResult[]                      │         Order[] │             │
│    │           │                          │                  │             │
│    │     Symptom[]                        │ (rsid link)     │             │
│    │           │                          │                  │             │
│    │     Treatment[]                      │                  │             │
│    │                                      │                  │             │
│    └─────► Practitioner ──────────────────│──────────────────┘             │
│                                           │                               │
│  ─────────────────────────────────────────│──────────────────────────── │
│                                           ▼                               │
│  RuVector Graph (Knowledge Domain - @ruvector/graph-node)                │
│  ═══════════════════════════════════════════════════════                │
│                                                                          │
│  Gene ◄────HAS_VARIANT───► SNP ◄────AFFECTS────► Nutrient               │
│    │                         │                       │                   │
│    │                         │                       │                   │
│  INVOLVED_IN            ASSOCIATED_WITH         CONTAINS                │
│    │                         │                       │                   │
│    ▼                         ▼                       ▼                   │
│  Pathway ◄─────────────► Condition ◄────────────► Food                  │
│    │                         │                                           │
│    │                         │                                           │
│  REQUIRES                  TREATS                                        │
│    │                         │                                           │
│    ▼                         ▼                                           │
│  Nutrient ◄───────────────► Herb ◄────INGREDIENT_OF────► Formula        │
│                              │                                           │
│                              │                                           │
│                         HAS_RESEARCH                                     │
│                              │                                           │
│                              ▼                                           │
│                         Publication                                      │
│                                                                          │
│  ╔═══════════════════════════════════════════════════════════════════╗  │
│  ║  HYPEREDGES (unique to RuVector Graph):                           ║  │
│  ║  • METHYLATION_CYCLE_PARTICIPANTS: [MTHFR, MTR, MTRR, SHMT1]     ║  │
│  ║  • FORMULA_COMPOSITION: [huang_qi, dang_shen, bai_zhu, gan_cao]  ║  │
│  ╚═══════════════════════════════════════════════════════════════════╝  │
│                                                                          │
│  ────────────────────────────────────────────────────────────────────── │
│                                                                          │
│  RuVector PostgreSQL (Vector Domain - Port 5432)                         │
│  ═══════════════════════════════════════════════                        │
│                                                                          │
│  embeddings ───── patterns ───── hyperbolic_embeddings                   │
│       │               │                    │                             │
│   HNSW Index    ReasoningBank        Poincare Ball                       │
│  (150x faster)   (AI Memory)      (Hierarchical Data)                    │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## Data Validation Rules

### User Domain

| Entity | Field | Validation |
|--------|-------|------------|
| User | email | Valid email format, unique |
| User | date_of_birth | Not in future, reasonable range |
| GeneticProfile | source | Enum: 23andme, ancestry, vcf, manual |
| UserSNP | rsid | Format: rs[0-9]+ |
| UserSNP | genotype | 1-2 characters, valid alleles |
| LabResult | value | Numeric, within reasonable bounds |
| Symptom | severity | Integer 1-10 |
| Treatment | effectiveness_rating | Integer 1-5 |

### Knowledge Domain

| Entity | Field | Validation |
|--------|-------|------------|
| Gene | symbol | HGNC approved symbol |
| SNP | rsid | Valid dbSNP identifier |
| Pathway | reactome_id | Valid Reactome accession |
| Herb | modality | Enum: tcm, ayurveda, kampo, western |
| Publication | pmid | Valid PubMed ID |

---

## Dependencies

| Document | Relationship |
|----------|--------------|
| [44-ARCHITECTURE](./44-ARCHITECTURE.md) | Technical architecture |
| [43-DATA-SOURCES](./43-DATA-SOURCES.md) | External data sources |
| [74-COMPLIANCE](../70-operations/74-COMPLIANCE.md) | Data compliance |

---

## Open Questions

- [ ] Define data retention policy for genetic data
- [ ] Establish anonymization strategy for research
- [ ] Determine backup and disaster recovery procedures
- [ ] Plan for GDPR right-to-deletion implementation

---

## Change Log

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | January 2026 | Engineering | Complete data model specification |
| 1.1 | January 2026 | Claude | Replace Neo4j with RuVector Graph; add hyperedge support |
