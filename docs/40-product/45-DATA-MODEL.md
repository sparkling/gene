# Core Data Model

**Document ID:** 45-DATA-MODEL
**Status:** Final
**Owner:** Engineering
**Last Updated:** January 2026
**Version:** 1.0

---

## TL;DR

Hybrid data model using PostgreSQL for user data/transactions, Neo4j for knowledge graph relationships. Core entities: User, GeneticProfile, SNP, Gene, Pathway, Nutrient, Herb, Condition, Practitioner. The knowledge graph connects THREE WORLDS: genes/SNPs (genetics) ↔ herbs/formulas (traditional medicine) ↔ nutrients/foods (nutrition). Designed for HIPAA compliance with field-level encryption on sensitive genetic data.

---

## Key Decisions

| Decision | Choice | Rationale | Date |
|----------|--------|-----------|------|
| User data storage | PostgreSQL | ACID, relational integrity | Jan 2026 |
| Knowledge relationships | Neo4j | Graph traversal performance | Jan 2026 |
| Genetic data encryption | Field-level AES-256 | HIPAA requirement | Jan 2026 |
| ID strategy | UUIDs | Distribution-friendly, privacy | Jan 2026 |

---

## Data Model Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                      USER DOMAIN (PostgreSQL)                    │
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
                              │ Links via external IDs
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    KNOWLEDGE DOMAIN (Neo4j)                      │
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

## Neo4j Schema (Knowledge Domain)

### Node Types

#### Gene Node

```cypher
CREATE (g:Gene {
    id: 'HGNC:7436',
    symbol: 'MTHFR',
    name: 'Methylenetetrahydrofolate reductase',
    chromosome: '1',
    description: 'Enzyme in folate metabolism',
    ensembl_id: 'ENSG00000177000',
    uniprot_id: 'P42898'
})
```

#### SNP Node

```cypher
CREATE (s:SNP {
    id: 'rs1801133',
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
})
```

#### Pathway Node

```cypher
CREATE (p:Pathway {
    id: 'REACT_R-HSA-196741',
    name: 'Methylation Cycle',
    description: 'One-carbon metabolism and methyl group transfer',
    category: 'metabolism',
    reactome_id: 'R-HSA-196741',
    kegg_id: 'hsa00670'
})
```

#### Nutrient Node

```cypher
CREATE (n:Nutrient {
    id: 'NUT_methylfolate',
    name: 'Methylfolate',
    alternate_names: ['5-MTHF', 'L-methylfolate', 'Metafolin'],
    category: 'vitamin',
    subcategory: 'B_vitamin',
    usda_id: '1177',
    recommended_daily: '400-800mcg'
})
```

#### Herb Node

```cypher
CREATE (h:Herb {
    id: 'TCM_huangqi',
    name: 'Huang Qi',
    latin_name: 'Astragalus membranaceus',
    english_name: 'Astragalus',
    modality: 'tcm',
    tcmsp_id: 'MOL000239',
    properties: ['qi_tonic', 'immune_support'],
    meridians: ['spleen', 'lung'],
    temperature: 'warm',
    taste: 'sweet'
})
```

#### Condition Node

```cypher
CREATE (c:Condition {
    id: 'COND_me_cfs',
    name: 'ME/CFS',
    full_name: 'Myalgic Encephalomyelitis / Chronic Fatigue Syndrome',
    icd10_code: 'G93.32',
    mesh_id: 'D015673',
    category: 'neurological'
})
```

#### Publication Node

```cypher
CREATE (pub:Publication {
    id: 'PMID_12345678',
    pmid: '12345678',
    title: 'MTHFR polymorphisms and folate status',
    authors: ['Smith J', 'Jones A'],
    journal: 'Nature Genetics',
    year: 2023,
    doi: '10.1038/ng.1234',
    evidence_level: 'clinical_trial'
})
```

#### Formula Node (TCM/Kampo)

```cypher
CREATE (f:Formula {
    id: 'FORMULA_buzhongyiqi',
    name: 'Bu Zhong Yi Qi Tang',
    english_name: 'Tonify the Middle and Augment the Qi Decoction',
    modality: 'tcm',
    category: 'qi_tonifying',
    indication: 'spleen_qi_deficiency'
})
```

### Relationship Types

```cypher
// Gene ←→ SNP
(gene:Gene)-[:HAS_VARIANT]->(snp:SNP)

// Gene ←→ Pathway
(gene:Gene)-[:INVOLVED_IN {role: 'enzyme'}]->(pathway:Pathway)

// SNP ←→ Nutrient
(snp:SNP)-[:AFFECTS {effect: 'reduces_activity', magnitude: 'moderate'}]->(nutrient:Nutrient)

// Pathway ←→ Nutrient
(pathway:Pathway)-[:REQUIRES_COFACTOR]->(nutrient:Nutrient)

// Herb ←→ Condition
(herb:Herb)-[:TREATS {evidence_level: 'clinical', mechanism: 'anti_inflammatory'}]->(condition:Condition)

// Herb ←→ Formula
(herb:Herb)-[:INGREDIENT_OF {role: 'chief'}]->(formula:Formula)

// SNP ←→ Condition
(snp:SNP)-[:ASSOCIATED_WITH {odds_ratio: 1.5, p_value: 0.001}]->(condition:Condition)

// Publication ←→ Various
(snp:SNP)-[:HAS_RESEARCH]->(pub:Publication)
(herb:Herb)-[:HAS_RESEARCH]->(pub:Publication)

// Food ←→ Nutrient
(food:Food)-[:CONTAINS {amount: 50, unit: 'mg', per: '100g'}]->(nutrient:Nutrient)
```

### Sample Queries

#### Find Nutrients Affected by User's SNPs

```cypher
MATCH (snp:SNP)-[a:AFFECTS]->(nutrient:Nutrient)
WHERE snp.rsid IN $user_rsids
RETURN snp.rsid, snp.common_name,
       nutrient.name, a.effect, a.magnitude
ORDER BY a.magnitude DESC
```

#### Find TCM Herbs for a Condition via Pathways

```cypher
MATCH (condition:Condition {name: $condition_name})
      <-[:TREATS]-(herb:Herb {modality: 'tcm'})
      -[:HAS_RESEARCH]->(pub:Publication)
WHERE pub.evidence_level IN ['clinical_trial', 'meta_analysis']
RETURN herb.name, herb.properties, count(pub) as research_count
ORDER BY research_count DESC
LIMIT 10
```

#### Find Pathway Genes with User's Risk Variants

```cypher
MATCH (pathway:Pathway {name: $pathway_name})
      <-[:INVOLVED_IN]-(gene:Gene)
      -[:HAS_VARIANT]->(snp:SNP)
WHERE snp.rsid IN $user_risk_snps
RETURN pathway.name, gene.symbol, snp.rsid, snp.clinical_significance
```

---

## Entity Relationship Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                    COMPLETE DATA MODEL                                   │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  PostgreSQL (User Domain)                                                │
│  ═══════════════════════                                                │
│                                                                          │
│  User ─────────┬────────────────────────────────────────┐               │
│    │           │                                         │               │
│    │     GeneticProfile                            Subscription         │
│    │           │                                         │               │
│    │       UserSNP[]────────────────────────────────────┼─┐             │
│    │           │                                         │ │             │
│    │     LabResult[]                                Order[] │             │
│    │           │                                             │             │
│    │     Symptom[]                                           │             │
│    │           │                                             │             │
│    │     Treatment[]                                         │             │
│    │                                                         │             │
│    └─────► Practitioner ─────────────────────────────────────┘             │
│                                                                          │
│  ────────────────────────────────────────────────────────────────────── │
│                                                                          │
│  Neo4j (Knowledge Domain)                                                │
│  ═══════════════════════                                                │
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
