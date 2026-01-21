/**
 * Cross-Database Query Service
 *
 * Bridges PostgreSQL (user data) with RuVector Graph Database (knowledge graph)
 * Enables queries that span both systems:
 * - User genetic data (PostgreSQL) -> Related knowledge (Graph)
 * - User symptoms/conditions -> Treatment recommendations (Graph)
 * - User lab results -> Nutrient/pathway analysis (Graph)
 *
 * Architecture:
 * - PostgreSQL: User profiles, genetic data, lab results, treatments
 * - RuVector Graph: THREE WORLDS knowledge (genetics, traditional medicine, nutrition)
 */

import { Pool, PoolClient, QueryResult } from 'pg';

// ============================================
// Types
// ============================================

interface PostgresConfig {
  host: string;
  port: number;
  database: string;
  user: string;
  password: string;
  max?: number;
  idleTimeoutMillis?: number;
}

interface GraphServiceConfig {
  baseUrl: string;
  timeout?: number;
}

interface UserSNP {
  id: string;
  rsid: string;
  genotype: string;
  riskLevel: string;
  chromosome?: string;
  position?: number;
}

interface UserLabResult {
  id: string;
  testName: string;
  value: number;
  unit: string;
  status: 'normal' | 'low' | 'high' | 'critical';
  category: string;
  testDate: Date;
}

interface UserTreatment {
  id: string;
  name: string;
  type: string;
  modality: string;
  dosage: string;
  effectivenessRating?: number;
}

interface GraphNode {
  id: string;
  labels: string[];
  properties: Record<string, string>;
}

interface GraphEdge {
  id: string;
  from: string;
  to: string;
  edgeType: string;
  properties: Record<string, string>;
}

interface GraphQueryResult {
  nodes: GraphNode[];
  edges: GraphEdge[];
  stats?: {
    totalNodes: number;
    totalEdges: number;
    avgDegree: number;
  };
}

interface NutrientRecommendation {
  nutrientId: string;
  nutrientName: string;
  reason: string;
  affectedGene: string;
  snpId: string;
  confidence: number;
  recommendedDose?: string;
}

interface HerbRecommendation {
  herbId: string;
  herbName: string;
  latinName: string;
  modality: string;
  targetCondition: string;
  confidence: number;
  pathwayTargets: string[];
}

interface UserHealthProfile {
  userId: string;
  riskSNPs: UserSNP[];
  abnormalLabs: UserLabResult[];
  currentTreatments: UserTreatment[];
  conditions: string[];
}

interface IntegratedRecommendation {
  userId: string;
  nutrients: NutrientRecommendation[];
  herbs: HerbRecommendation[];
  pathwayRisks: string[];
  generatedAt: Date;
}

// ============================================
// PostgreSQL Client
// ============================================

class PostgresClient {
  private pool: Pool;

  constructor(config: PostgresConfig) {
    this.pool = new Pool({
      host: config.host,
      port: config.port,
      database: config.database,
      user: config.user,
      password: config.password,
      max: config.max || 20,
      idleTimeoutMillis: config.idleTimeoutMillis || 30000,
    });

    this.pool.on('error', (err) => {
      console.error('[PostgresClient] Pool error:', err);
    });
  }

  async query<T = any>(sql: string, params?: any[]): Promise<QueryResult<T>> {
    const client = await this.pool.connect();
    try {
      return await client.query(sql, params);
    } finally {
      client.release();
    }
  }

  async getClient(): Promise<PoolClient> {
    return this.pool.connect();
  }

  async close(): Promise<void> {
    await this.pool.end();
  }

  // ==========================================
  // User Data Queries
  // ==========================================

  async getUserSNPs(userId: string): Promise<UserSNP[]> {
    const result = await this.query<UserSNP>(`
      SELECT
        us.id,
        us.rsid,
        us.genotype,
        us.risk_level as "riskLevel",
        us.chromosome,
        us.position
      FROM user_snps us
      JOIN genetic_profiles gp ON us.genetic_profile_id = gp.id
      WHERE gp.user_id = $1
        AND gp.status = 'complete'
      ORDER BY us.risk_level DESC, us.rsid
    `, [userId]);

    return result.rows;
  }

  async getUserRiskSNPs(userId: string): Promise<UserSNP[]> {
    const result = await this.query<UserSNP>(`
      SELECT
        us.id,
        us.rsid,
        us.genotype,
        us.risk_level as "riskLevel",
        us.chromosome,
        us.position
      FROM user_snps us
      JOIN genetic_profiles gp ON us.genetic_profile_id = gp.id
      WHERE gp.user_id = $1
        AND gp.status = 'complete'
        AND us.risk_level IN ('heterozygous', 'homozygous_risk')
      ORDER BY
        CASE us.risk_level
          WHEN 'homozygous_risk' THEN 1
          WHEN 'heterozygous' THEN 2
          ELSE 3
        END,
        us.rsid
    `, [userId]);

    return result.rows;
  }

  async getUserLabResults(userId: string, daysBack: number = 90): Promise<UserLabResult[]> {
    const result = await this.query<UserLabResult>(`
      SELECT
        id,
        test_name as "testName",
        value,
        unit,
        status,
        category,
        test_date as "testDate"
      FROM lab_results
      WHERE user_id = $1
        AND test_date >= NOW() - INTERVAL '${daysBack} days'
      ORDER BY test_date DESC, test_name
    `, [userId]);

    return result.rows;
  }

  async getUserAbnormalLabs(userId: string): Promise<UserLabResult[]> {
    const result = await this.query<UserLabResult>(`
      SELECT
        id,
        test_name as "testName",
        value,
        unit,
        status,
        category,
        test_date as "testDate"
      FROM lab_results
      WHERE user_id = $1
        AND status IN ('low', 'high', 'critical')
        AND test_date >= NOW() - INTERVAL '90 days'
      ORDER BY
        CASE status
          WHEN 'critical' THEN 1
          WHEN 'high' THEN 2
          WHEN 'low' THEN 3
          ELSE 4
        END,
        test_date DESC
    `, [userId]);

    return result.rows;
  }

  async getUserTreatments(userId: string, activeOnly: boolean = true): Promise<UserTreatment[]> {
    let sql = `
      SELECT
        id,
        name,
        type,
        modality,
        dosage,
        effectiveness_rating as "effectivenessRating"
      FROM treatments
      WHERE user_id = $1
    `;

    if (activeOnly) {
      sql += ` AND (ended_at IS NULL OR ended_at > NOW())`;
    }

    sql += ` ORDER BY started_at DESC`;

    const result = await this.query<UserTreatment>(sql, [userId]);
    return result.rows;
  }

  async getUserConditions(userId: string): Promise<string[]> {
    // Get conditions from symptoms that have been mapped to conditions
    const result = await this.query<{ name: string }>(`
      SELECT DISTINCT name
      FROM symptoms
      WHERE user_id = $1
        AND (ended_at IS NULL OR ended_at > NOW())
        AND severity >= 5
      ORDER BY name
    `, [userId]);

    return result.rows.map(r => r.name);
  }
}

// ============================================
// Graph Service Client
// ============================================

class GraphServiceClient {
  private baseUrl: string;
  private timeout: number;

  constructor(config: GraphServiceConfig) {
    this.baseUrl = config.baseUrl.replace(/\/$/, '');
    this.timeout = config.timeout || 30000;
  }

  private async fetch<T>(path: string, options: RequestInit = {}): Promise<T> {
    const controller = new AbortController();
    const timeoutId = setTimeout(() => controller.abort(), this.timeout);

    try {
      const response = await fetch(`${this.baseUrl}${path}`, {
        ...options,
        signal: controller.signal,
        headers: {
          'Content-Type': 'application/json',
          ...options.headers,
        },
      });

      if (!response.ok) {
        const error = await response.text();
        throw new Error(`Graph service error: ${response.status} - ${error}`);
      }

      return response.json();
    } finally {
      clearTimeout(timeoutId);
    }
  }

  async query(cypher: string): Promise<GraphQueryResult> {
    return this.fetch('/query', {
      method: 'POST',
      body: JSON.stringify({ cypher }),
    });
  }

  async getNeighbors(nodeId: string, k: number = 2): Promise<string[]> {
    const result = await this.fetch<{ neighbors: string[] }>(`/neighbors/${nodeId}?k=${k}`);
    return result.neighbors;
  }

  async searchHyperedges(embedding: number[], k: number = 10): Promise<any[]> {
    const result = await this.fetch<{ results: any[] }>('/hyperedges/search', {
      method: 'POST',
      body: JSON.stringify({ embedding, k }),
    });
    return result.results;
  }

  // ==========================================
  // THREE WORLDS Queries
  // ==========================================

  async getNutrientsAffectedBySNP(rsid: string): Promise<GraphNode[]> {
    const result = await this.query(`
      MATCH (s:SNP)-[a:AFFECTS]->(n:Nutrient)
      WHERE s.rsid = '${rsid}'
      RETURN n, a
    `);
    return result.nodes.filter(n => n.labels.includes('Nutrient'));
  }

  async getPathwaysByGene(geneSymbol: string): Promise<GraphNode[]> {
    const result = await this.query(`
      MATCH (g:Gene)-[:INVOLVED_IN]->(p:Pathway)
      WHERE g.symbol = '${geneSymbol}'
      RETURN p
    `);
    return result.nodes.filter(n => n.labels.includes('Pathway'));
  }

  async getHerbsForCondition(conditionName: string): Promise<GraphNode[]> {
    const result = await this.query(`
      MATCH (h:Herb)-[:TREATS]->(c:Condition)
      WHERE c.name CONTAINS '${conditionName}'
      RETURN h
    `);
    return result.nodes.filter(n => n.labels.includes('Herb'));
  }

  async getGeneFromSNP(rsid: string): Promise<GraphNode | null> {
    const result = await this.query(`
      MATCH (g:Gene)-[:HAS_VARIANT]->(s:SNP)
      WHERE s.rsid = '${rsid}'
      RETURN g
    `);
    const genes = result.nodes.filter(n => n.labels.includes('Gene'));
    return genes[0] || null;
  }

  async getCofactorsForPathway(pathwayId: string): Promise<GraphNode[]> {
    const result = await this.query(`
      MATCH (p:Pathway)-[:REQUIRES_COFACTOR]->(n:Nutrient)
      WHERE p.id = '${pathwayId}'
      RETURN n
    `);
    return result.nodes.filter(n => n.labels.includes('Nutrient'));
  }

  async getConditionsFromSNP(rsid: string): Promise<GraphNode[]> {
    const result = await this.query(`
      MATCH (s:SNP)-[:ASSOCIATED_WITH]->(c:Condition)
      WHERE s.rsid = '${rsid}'
      RETURN c
    `);
    return result.nodes.filter(n => n.labels.includes('Condition'));
  }
}

// ============================================
// Cross-Database Query Service
// ============================================

export class CrossDatabaseService {
  private postgres: PostgresClient;
  private graph: GraphServiceClient;

  constructor(postgresConfig: PostgresConfig, graphConfig: GraphServiceConfig) {
    this.postgres = new PostgresClient(postgresConfig);
    this.graph = new GraphServiceClient(graphConfig);
  }

  async close(): Promise<void> {
    await this.postgres.close();
  }

  // ==========================================
  // User Health Profile
  // ==========================================

  async getUserHealthProfile(userId: string): Promise<UserHealthProfile> {
    const [riskSNPs, abnormalLabs, currentTreatments, conditions] = await Promise.all([
      this.postgres.getUserRiskSNPs(userId),
      this.postgres.getUserAbnormalLabs(userId),
      this.postgres.getUserTreatments(userId, true),
      this.postgres.getUserConditions(userId),
    ]);

    return {
      userId,
      riskSNPs,
      abnormalLabs,
      currentTreatments,
      conditions,
    };
  }

  // ==========================================
  // Integrated Recommendations
  // ==========================================

  async generateRecommendations(userId: string): Promise<IntegratedRecommendation> {
    const profile = await this.getUserHealthProfile(userId);

    // Generate nutrient recommendations based on SNPs
    const nutrientRecommendations: NutrientRecommendation[] = [];
    for (const snp of profile.riskSNPs) {
      const nutrients = await this.graph.getNutrientsAffectedBySNP(snp.rsid);
      const gene = await this.graph.getGeneFromSNP(snp.rsid);

      for (const nutrient of nutrients) {
        nutrientRecommendations.push({
          nutrientId: nutrient.id,
          nutrientName: nutrient.properties.name || nutrient.id,
          reason: `SNP ${snp.rsid} (${snp.genotype}) affects ${nutrient.properties.name} metabolism`,
          affectedGene: gene?.properties.symbol || 'Unknown',
          snpId: snp.rsid,
          confidence: snp.riskLevel === 'homozygous_risk' ? 0.9 : 0.7,
          recommendedDose: nutrient.properties.recommendedDaily,
        });
      }
    }

    // Generate herb recommendations based on conditions
    const herbRecommendations: HerbRecommendation[] = [];
    for (const conditionName of profile.conditions) {
      const herbs = await this.graph.getHerbsForCondition(conditionName);

      for (const herb of herbs) {
        herbRecommendations.push({
          herbId: herb.id,
          herbName: herb.properties.name || herb.id,
          latinName: herb.properties.latinName || '',
          modality: herb.properties.modality || 'unknown',
          targetCondition: conditionName,
          confidence: 0.7,
          pathwayTargets: [],
        });
      }
    }

    // Identify pathway risks from SNPs
    const pathwayRisks: string[] = [];
    for (const snp of profile.riskSNPs) {
      const gene = await this.graph.getGeneFromSNP(snp.rsid);
      if (gene) {
        const pathways = await this.graph.getPathwaysByGene(gene.properties.symbol);
        for (const pathway of pathways) {
          const riskDesc = `${pathway.properties.name} (affected by ${gene.properties.symbol} ${snp.rsid})`;
          if (!pathwayRisks.includes(riskDesc)) {
            pathwayRisks.push(riskDesc);
          }
        }
      }
    }

    return {
      userId,
      nutrients: nutrientRecommendations,
      herbs: herbRecommendations,
      pathwayRisks,
      generatedAt: new Date(),
    };
  }

  // ==========================================
  // SNP Analysis
  // ==========================================

  async analyzeSNPImpact(userId: string, rsid: string): Promise<{
    snp: UserSNP | null;
    gene: GraphNode | null;
    pathways: GraphNode[];
    affectedNutrients: GraphNode[];
    associatedConditions: GraphNode[];
    cofactors: GraphNode[];
  }> {
    // Get user's SNP data
    const userSNPs = await this.postgres.getUserSNPs(userId);
    const snp = userSNPs.find(s => s.rsid === rsid) || null;

    // Get graph data
    const [gene, affectedNutrients, associatedConditions] = await Promise.all([
      this.graph.getGeneFromSNP(rsid),
      this.graph.getNutrientsAffectedBySNP(rsid),
      this.graph.getConditionsFromSNP(rsid),
    ]);

    // Get pathways and cofactors if gene found
    let pathways: GraphNode[] = [];
    let cofactors: GraphNode[] = [];

    if (gene) {
      pathways = await this.graph.getPathwaysByGene(gene.properties.symbol);

      // Get cofactors for all pathways
      for (const pathway of pathways) {
        const pathwayCofactors = await this.graph.getCofactorsForPathway(pathway.id);
        cofactors.push(...pathwayCofactors);
      }

      // Deduplicate cofactors
      cofactors = cofactors.filter((c, i, arr) =>
        arr.findIndex(x => x.id === c.id) === i
      );
    }

    return {
      snp,
      gene,
      pathways,
      affectedNutrients,
      associatedConditions,
      cofactors,
    };
  }

  // ==========================================
  // Treatment Efficacy Analysis
  // ==========================================

  async analyzeTreatmentEfficacy(userId: string): Promise<{
    treatments: Array<{
      treatment: UserTreatment;
      targetedPathways: GraphNode[];
      targetedConditions: GraphNode[];
      potentialInteractions: string[];
    }>;
  }> {
    const treatments = await this.postgres.getUserTreatments(userId);
    const results: Array<{
      treatment: UserTreatment;
      targetedPathways: GraphNode[];
      targetedConditions: GraphNode[];
      potentialInteractions: string[];
    }> = [];

    for (const treatment of treatments) {
      let targetedPathways: GraphNode[] = [];
      let targetedConditions: GraphNode[] = [];
      const potentialInteractions: string[] = [];

      // For herbs, query the graph
      if (treatment.type === 'herb' || treatment.modality === 'tcm' || treatment.modality === 'ayurveda') {
        // Query herbs by name
        const herbResult = await this.graph.query(`
          MATCH (h:Herb)
          WHERE h.name CONTAINS '${treatment.name}' OR h.englishName CONTAINS '${treatment.name}'
          RETURN h
        `);

        if (herbResult.nodes.length > 0) {
          const herb = herbResult.nodes[0];

          // Get targeted conditions
          targetedConditions = await this.graph.getHerbsForCondition(treatment.name);

          // Get pathway targets
          const pathwayResult = await this.graph.query(`
            MATCH (h:Herb)-[:TARGETS]->(p:Pathway)
            WHERE h.id = '${herb.id}'
            RETURN p
          `);
          targetedPathways = pathwayResult.nodes.filter(n => n.labels.includes('Pathway'));
        }
      }

      // For supplements, check nutrient interactions
      if (treatment.type === 'supplement') {
        const nutrientResult = await this.graph.query(`
          MATCH (n:Nutrient)
          WHERE n.name CONTAINS '${treatment.name}'
          RETURN n
        `);

        if (nutrientResult.nodes.length > 0) {
          // Get pathways requiring this nutrient
          const pathwayResult = await this.graph.query(`
            MATCH (p:Pathway)-[:REQUIRES_COFACTOR]->(n:Nutrient)
            WHERE n.name CONTAINS '${treatment.name}'
            RETURN p
          `);
          targetedPathways = pathwayResult.nodes.filter(n => n.labels.includes('Pathway'));
        }
      }

      results.push({
        treatment,
        targetedPathways,
        targetedConditions,
        potentialInteractions,
      });
    }

    return { treatments: results };
  }

  // ==========================================
  // Lab Result Correlation
  // ==========================================

  async correlateLabsWithGenetics(userId: string): Promise<{
    correlations: Array<{
      labResult: UserLabResult;
      relatedSNPs: UserSNP[];
      relatedNutrients: GraphNode[];
      explanation: string;
    }>;
  }> {
    const [labResults, userSNPs] = await Promise.all([
      this.postgres.getUserAbnormalLabs(userId),
      this.postgres.getUserRiskSNPs(userId),
    ]);

    const correlations: Array<{
      labResult: UserLabResult;
      relatedSNPs: UserSNP[];
      relatedNutrients: GraphNode[];
      explanation: string;
    }> = [];

    // Map common lab tests to nutrients/pathways
    const labNutrientMap: Record<string, string[]> = {
      'homocysteine': ['Methylfolate', 'B12', 'B6'],
      'vitamin d': ['Vitamin D'],
      'b12': ['Vitamin B12'],
      'folate': ['Methylfolate', 'Folate'],
      'ferritin': ['Iron'],
      'magnesium': ['Magnesium'],
      'zinc': ['Zinc'],
    };

    for (const lab of labResults) {
      const labNameLower = lab.testName.toLowerCase();
      const relatedSNPs: UserSNP[] = [];
      const relatedNutrients: GraphNode[] = [];
      let explanation = '';

      // Find related nutrients
      const nutrientNames = Object.entries(labNutrientMap)
        .filter(([key]) => labNameLower.includes(key))
        .flatMap(([, nutrients]) => nutrients);

      for (const nutrientName of nutrientNames) {
        // Query graph for nutrient
        const nutrientResult = await this.graph.query(`
          MATCH (n:Nutrient)
          WHERE n.name CONTAINS '${nutrientName}'
          RETURN n
        `);
        relatedNutrients.push(...nutrientResult.nodes);

        // Find SNPs affecting this nutrient
        for (const snp of userSNPs) {
          const affectedNutrients = await this.graph.getNutrientsAffectedBySNP(snp.rsid);
          if (affectedNutrients.some(n => n.properties.name?.includes(nutrientName))) {
            if (!relatedSNPs.find(s => s.rsid === snp.rsid)) {
              relatedSNPs.push(snp);
            }
          }
        }
      }

      // Generate explanation
      if (relatedSNPs.length > 0) {
        const snpList = relatedSNPs.map(s => s.rsid).join(', ');
        explanation = `${lab.status === 'low' ? 'Low' : 'High'} ${lab.testName} may be related to genetic variants (${snpList}) affecting nutrient metabolism.`;
      } else {
        explanation = `${lab.testName} is ${lab.status}. Consider dietary and lifestyle factors.`;
      }

      correlations.push({
        labResult: lab,
        relatedSNPs,
        relatedNutrients,
        explanation,
      });
    }

    return { correlations };
  }
}

// ============================================
// Factory Function
// ============================================

export function createCrossDatabaseService(): CrossDatabaseService {
  const postgresConfig: PostgresConfig = {
    host: process.env.POSTGRES_HOST || 'localhost',
    port: parseInt(process.env.POSTGRES_PORT || '5433', 10),
    database: process.env.POSTGRES_DB || 'gene_user_db',
    user: process.env.POSTGRES_USER || 'gene_user',
    password: process.env.POSTGRES_PASSWORD || 'gene_secure_password',
  };

  const graphConfig: GraphServiceConfig = {
    baseUrl: process.env.GRAPH_SERVICE_URL || 'http://localhost:3100',
    timeout: parseInt(process.env.GRAPH_TIMEOUT || '30000', 10),
  };

  return new CrossDatabaseService(postgresConfig, graphConfig);
}

export default CrossDatabaseService;
