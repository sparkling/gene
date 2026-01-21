-- ============================================
-- USER DATABASE INITIALIZATION (PostgreSQL)
-- ============================================
--
-- This script initializes the PostgreSQL user database
-- with RuVector extension for vector search capabilities.
--
-- Schema includes:
-- - User profiles and authentication
-- - Genetic profiles and SNP data
-- - Lab results and symptoms
-- - Treatments and subscriptions
-- - Practitioner marketplace data
--
-- Designed for HIPAA compliance with field-level encryption

-- ============================================
-- EXTENSIONS
-- ============================================

-- RuVector for vector similarity search
CREATE EXTENSION IF NOT EXISTS ruvector VERSION '0.1.0';
CREATE EXTENSION IF NOT EXISTS pgcrypto;
CREATE EXTENSION IF NOT EXISTS "uuid-ossp";

-- ============================================
-- SCHEMA
-- ============================================

CREATE SCHEMA IF NOT EXISTS gene_user;
SET search_path TO gene_user, public;

-- ============================================
-- USER TABLES
-- ============================================

-- Users table
CREATE TABLE IF NOT EXISTS gene_user.users (
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

    -- Preferences (encrypted for sensitive data)
    preferences JSONB DEFAULT '{}',

    -- User embedding for personalized recommendations
    user_embedding ruvector(384),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    last_login_at TIMESTAMPTZ,
    deleted_at TIMESTAMPTZ -- soft delete
);

CREATE INDEX idx_users_email ON gene_user.users(email);
CREATE INDEX idx_users_role ON gene_user.users(role);
CREATE INDEX idx_users_embedding ON gene_user.users
    USING hnsw (user_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

-- Genetic Profiles
CREATE TABLE IF NOT EXISTS gene_user.genetic_profiles (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id) ON DELETE CASCADE,

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

    -- Encrypted raw data (HIPAA compliant)
    raw_data_encrypted BYTEA, -- AES-256 encrypted
    raw_data_iv BYTEA,

    -- Profile embedding for similarity matching
    profile_embedding ruvector(384),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_genetic_profiles_user ON gene_user.genetic_profiles(user_id);
CREATE INDEX idx_genetic_profiles_status ON gene_user.genetic_profiles(status);
CREATE INDEX idx_genetic_profiles_embedding ON gene_user.genetic_profiles
    USING hnsw (profile_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

-- User SNPs
CREATE TABLE IF NOT EXISTS gene_user.user_snps (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    genetic_profile_id UUID NOT NULL REFERENCES gene_user.genetic_profiles(id) ON DELETE CASCADE,

    -- SNP Data
    rsid VARCHAR(20) NOT NULL, -- e.g., 'rs1801133'
    chromosome VARCHAR(5),
    position BIGINT,
    genotype VARCHAR(10), -- e.g., 'CT', 'AA'

    -- Reference (links to Graph DB)
    graph_snp_id VARCHAR(100), -- Reference to RuVector graph node

    -- Interpretation
    risk_level VARCHAR(20), -- 'normal', 'heterozygous', 'homozygous_risk'
    clinical_significance VARCHAR(50),

    -- SNP embedding for similarity queries
    snp_embedding ruvector(384),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),

    CONSTRAINT unique_profile_rsid UNIQUE (genetic_profile_id, rsid)
);

CREATE INDEX idx_user_snps_profile ON gene_user.user_snps(genetic_profile_id);
CREATE INDEX idx_user_snps_rsid ON gene_user.user_snps(rsid);
CREATE INDEX idx_user_snps_risk ON gene_user.user_snps(risk_level);
CREATE INDEX idx_user_snps_embedding ON gene_user.user_snps
    USING hnsw (snp_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 100);

-- Lab Results
CREATE TABLE IF NOT EXISTS gene_user.lab_results (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id) ON DELETE CASCADE,

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

    -- Reference to graph biomarkers
    graph_biomarker_id VARCHAR(100),

    -- Embedding for semantic search
    result_embedding ruvector(384),

    -- Metadata
    test_date DATE,
    lab_name VARCHAR(255),
    notes TEXT,

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_lab_results_user ON gene_user.lab_results(user_id);
CREATE INDEX idx_lab_results_category ON gene_user.lab_results(category);
CREATE INDEX idx_lab_results_date ON gene_user.lab_results(test_date);
CREATE INDEX idx_lab_results_status ON gene_user.lab_results(status);
CREATE INDEX idx_lab_results_embedding ON gene_user.lab_results
    USING hnsw (result_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

-- Symptoms
CREATE TABLE IF NOT EXISTS gene_user.symptoms (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id) ON DELETE CASCADE,

    -- Symptom
    name VARCHAR(255) NOT NULL,
    category VARCHAR(50), -- 'physical', 'mental', 'digestive', 'neurological'

    -- Reference to graph conditions
    graph_condition_id VARCHAR(100),

    -- Severity & Tracking
    severity INTEGER CHECK (severity BETWEEN 1 AND 10),
    frequency VARCHAR(50), -- 'constant', 'daily', 'weekly', 'occasional'

    -- Timing
    started_at DATE,
    ended_at DATE,

    -- Notes
    notes TEXT,
    triggers TEXT[],

    -- Embedding for symptom clustering
    symptom_embedding ruvector(384),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_symptoms_user ON gene_user.symptoms(user_id);
CREATE INDEX idx_symptoms_category ON gene_user.symptoms(category);
CREATE INDEX idx_symptoms_severity ON gene_user.symptoms(severity);
CREATE INDEX idx_symptoms_embedding ON gene_user.symptoms
    USING hnsw (symptom_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

-- Treatments
CREATE TABLE IF NOT EXISTS gene_user.treatments (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id) ON DELETE CASCADE,

    -- Treatment
    name VARCHAR(255) NOT NULL,
    type VARCHAR(50), -- 'supplement', 'medication', 'herb', 'lifestyle'
    modality VARCHAR(50), -- 'western', 'tcm', 'ayurveda', 'kampo'

    -- Reference to graph entities
    graph_herb_id VARCHAR(100),
    graph_nutrient_id VARCHAR(100),
    graph_formula_id VARCHAR(100),

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

    -- Embedding for treatment similarity
    treatment_embedding ruvector(384),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_treatments_user ON gene_user.treatments(user_id);
CREATE INDEX idx_treatments_type ON gene_user.treatments(type);
CREATE INDEX idx_treatments_modality ON gene_user.treatments(modality);
CREATE INDEX idx_treatments_active ON gene_user.treatments(ended_at) WHERE ended_at IS NULL;
CREATE INDEX idx_treatments_embedding ON gene_user.treatments
    USING hnsw (treatment_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

-- Subscriptions
CREATE TABLE IF NOT EXISTS gene_user.subscriptions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id) ON DELETE CASCADE,

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

CREATE INDEX idx_subscriptions_user ON gene_user.subscriptions(user_id);
CREATE INDEX idx_subscriptions_status ON gene_user.subscriptions(status);
CREATE INDEX idx_subscriptions_stripe ON gene_user.subscriptions(stripe_subscription_id);

-- Orders (Marketplace)
CREATE TABLE IF NOT EXISTS gene_user.orders (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id),

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

CREATE INDEX idx_orders_user ON gene_user.orders(user_id);
CREATE INDEX idx_orders_status ON gene_user.orders(status);
CREATE INDEX idx_orders_created ON gene_user.orders(created_at);

-- Practitioners
CREATE TABLE IF NOT EXISTS gene_user.practitioners (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id) ON DELETE CASCADE,

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

    -- Practitioner embedding for matching
    practitioner_embedding ruvector(384),

    -- Status
    status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'approved', 'suspended'

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),

    CONSTRAINT unique_practitioner_user UNIQUE (user_id)
);

CREATE INDEX idx_practitioners_status ON gene_user.practitioners(status);
CREATE INDEX idx_practitioners_specialties ON gene_user.practitioners USING GIN(specialties);
CREATE INDEX idx_practitioners_modalities ON gene_user.practitioners USING GIN(modalities);
CREATE INDEX idx_practitioners_embedding ON gene_user.practitioners
    USING hnsw (practitioner_embedding ruvector_cosine_ops)
    WITH (m = 16, ef_construction = 64);

-- ============================================
-- AUDIT & COMPLIANCE TABLES
-- ============================================

-- Audit log for HIPAA compliance
CREATE TABLE IF NOT EXISTS gene_user.audit_log (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES gene_user.users(id),
    action VARCHAR(50) NOT NULL, -- 'view', 'create', 'update', 'delete', 'export'
    resource_type VARCHAR(50) NOT NULL, -- 'genetic_profile', 'lab_result', etc.
    resource_id UUID,
    ip_address INET,
    user_agent TEXT,
    details JSONB,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX idx_audit_log_user ON gene_user.audit_log(user_id);
CREATE INDEX idx_audit_log_action ON gene_user.audit_log(action);
CREATE INDEX idx_audit_log_created ON gene_user.audit_log(created_at);

-- Data export requests (GDPR compliance)
CREATE TABLE IF NOT EXISTS gene_user.data_export_requests (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id),
    status VARCHAR(20) DEFAULT 'pending', -- 'pending', 'processing', 'complete', 'failed'
    file_url VARCHAR(500),
    expires_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    completed_at TIMESTAMPTZ
);

CREATE INDEX idx_data_exports_user ON gene_user.data_export_requests(user_id);
CREATE INDEX idx_data_exports_status ON gene_user.data_export_requests(status);

-- ============================================
-- GRAPH LINKAGE TABLE
-- ============================================

-- Maps user entities to RuVector graph nodes
CREATE TABLE IF NOT EXISTS gene_user.graph_links (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID NOT NULL REFERENCES gene_user.users(id),
    entity_type VARCHAR(50) NOT NULL, -- 'snp', 'treatment', 'symptom', 'lab_result'
    entity_id UUID NOT NULL,
    graph_node_id VARCHAR(255) NOT NULL, -- ID in RuVector graph
    graph_node_type VARCHAR(50) NOT NULL, -- 'Gene', 'SNP', 'Herb', 'Nutrient', etc.
    link_confidence FLOAT DEFAULT 1.0,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    UNIQUE(entity_type, entity_id, graph_node_id)
);

CREATE INDEX idx_graph_links_user ON gene_user.graph_links(user_id);
CREATE INDEX idx_graph_links_entity ON gene_user.graph_links(entity_type, entity_id);
CREATE INDEX idx_graph_links_graph ON gene_user.graph_links(graph_node_id);

-- ============================================
-- SEARCH FUNCTIONS
-- ============================================

-- Find similar users by embedding
CREATE OR REPLACE FUNCTION gene_user.find_similar_users(
    query_embedding ruvector(384),
    limit_count INT DEFAULT 10,
    min_similarity FLOAT DEFAULT 0.5
)
RETURNS TABLE (
    user_id UUID,
    email VARCHAR(255),
    similarity FLOAT
) AS $$
BEGIN
    RETURN QUERY
    SELECT
        u.id,
        u.email,
        (1 - (u.user_embedding <=> query_embedding))::FLOAT AS similarity
    FROM gene_user.users u
    WHERE u.user_embedding IS NOT NULL
      AND u.deleted_at IS NULL
      AND (1 - (u.user_embedding <=> query_embedding)) >= min_similarity
    ORDER BY u.user_embedding <=> query_embedding
    LIMIT limit_count;
END;
$$ LANGUAGE plpgsql STABLE;

-- Search SNPs by embedding similarity
CREATE OR REPLACE FUNCTION gene_user.search_snps_by_embedding(
    query_embedding ruvector(384),
    user_id_filter UUID DEFAULT NULL,
    limit_count INT DEFAULT 10
)
RETURNS TABLE (
    snp_id UUID,
    rsid VARCHAR(20),
    genotype VARCHAR(10),
    risk_level VARCHAR(20),
    similarity FLOAT
) AS $$
BEGIN
    RETURN QUERY
    SELECT
        us.id,
        us.rsid,
        us.genotype,
        us.risk_level,
        (1 - (us.snp_embedding <=> query_embedding))::FLOAT AS similarity
    FROM gene_user.user_snps us
    JOIN gene_user.genetic_profiles gp ON us.genetic_profile_id = gp.id
    WHERE us.snp_embedding IS NOT NULL
      AND (user_id_filter IS NULL OR gp.user_id = user_id_filter)
    ORDER BY us.snp_embedding <=> query_embedding
    LIMIT limit_count;
END;
$$ LANGUAGE plpgsql STABLE;

-- Find practitioners by specialty embedding
CREATE OR REPLACE FUNCTION gene_user.find_practitioners(
    query_embedding ruvector(384),
    modality_filter VARCHAR(50) DEFAULT NULL,
    limit_count INT DEFAULT 10
)
RETURNS TABLE (
    practitioner_id UUID,
    display_name VARCHAR(255),
    specialties TEXT[],
    modalities TEXT[],
    similarity FLOAT
) AS $$
BEGIN
    RETURN QUERY
    SELECT
        p.id,
        p.display_name,
        p.specialties,
        p.modalities,
        (1 - (p.practitioner_embedding <=> query_embedding))::FLOAT AS similarity
    FROM gene_user.practitioners p
    WHERE p.practitioner_embedding IS NOT NULL
      AND p.status = 'approved'
      AND p.accepting_new_clients = TRUE
      AND (modality_filter IS NULL OR modality_filter = ANY(p.modalities))
    ORDER BY p.practitioner_embedding <=> query_embedding
    LIMIT limit_count;
END;
$$ LANGUAGE plpgsql STABLE;

-- ============================================
-- COMPLETION
-- ============================================

DO $$
BEGIN
    RAISE NOTICE '';
    RAISE NOTICE '============================================';
    RAISE NOTICE 'User Database Initialization Complete!';
    RAISE NOTICE '============================================';
    RAISE NOTICE '';
    RAISE NOTICE 'Schema: gene_user';
    RAISE NOTICE '';
    RAISE NOTICE 'Tables:';
    RAISE NOTICE '  - users (with RuVector embedding)';
    RAISE NOTICE '  - genetic_profiles';
    RAISE NOTICE '  - user_snps (with RuVector embedding)';
    RAISE NOTICE '  - lab_results (with RuVector embedding)';
    RAISE NOTICE '  - symptoms (with RuVector embedding)';
    RAISE NOTICE '  - treatments (with RuVector embedding)';
    RAISE NOTICE '  - subscriptions';
    RAISE NOTICE '  - orders';
    RAISE NOTICE '  - practitioners (with RuVector embedding)';
    RAISE NOTICE '  - audit_log';
    RAISE NOTICE '  - data_export_requests';
    RAISE NOTICE '  - graph_links';
    RAISE NOTICE '';
    RAISE NOTICE 'HNSW Indices: 7 vector indices for similarity search';
    RAISE NOTICE 'Graph Links: Maps user data to RuVector graph nodes';
    RAISE NOTICE '';
END $$;
