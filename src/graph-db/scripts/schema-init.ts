/**
 * THREE WORLDS Schema Initialization
 *
 * Initializes the RuVector graph database with schema for:
 * - World 1: Genetics (Gene, SNP, Pathway, Condition, Publication)
 * - World 2: Traditional Medicine (Herb, Formula, Modality, TCMPattern, Meridian)
 * - World 3: Nutrition (Nutrient, Food, Biomarker, Supplement, FoodGroup)
 *
 * Uses @ruvector/graph-node native bindings
 */

const { GraphDatabase } = require('@ruvector/graph-node');

// ============================================
// Configuration
// ============================================

const config = {
  storagePath: process.env.GRAPH_STORAGE_PATH || './data/knowledge.db',
  dimensions: parseInt(process.env.GRAPH_DIMENSIONS || '384', 10),
  distanceMetric: 'Cosine' as const,
};

// ============================================
// Helper Functions
// ============================================

/**
 * Generate a zero embedding for schema nodes
 * Real embeddings will be populated during data import
 */
function zeroEmbedding(dim: number = 384): Float32Array {
  return new Float32Array(dim);
}

/**
 * Generate a random embedding for testing
 */
function randomEmbedding(dim: number = 384): Float32Array {
  const embedding = new Float32Array(dim);
  for (let i = 0; i < dim; i++) {
    embedding[i] = (Math.random() - 0.5) * 2;
  }
  // Normalize
  let norm = 0;
  for (let i = 0; i < dim; i++) {
    norm += embedding[i] * embedding[i];
  }
  norm = Math.sqrt(norm);
  for (let i = 0; i < dim; i++) {
    embedding[i] /= norm;
  }
  return embedding;
}

// ============================================
// Schema Definition: World 1 - Genetics
// ============================================

interface GeneNode {
  id: string;
  symbol: string;
  name: string;
  chromosome?: string;
  description?: string;
  ensemblId?: string;
  uniprotId?: string;
  hgncId?: string;
}

interface SNPNode {
  id: string;
  rsid: string;
  geneSymbol?: string;
  chromosome?: string;
  position?: number;
  refAllele?: string;
  altAllele?: string;
  minorAllele?: string;
  maf?: number;
  clinicalSignificance?: string;
  commonName?: string;
}

interface PathwayNode {
  id: string;
  name: string;
  description?: string;
  category?: string;
  reactomeId?: string;
  keggId?: string;
}

interface ConditionNode {
  id: string;
  name: string;
  fullName?: string;
  icd10Code?: string;
  meshId?: string;
  category?: string;
  mondoId?: string;
  orphanetId?: string;
}

interface PublicationNode {
  id: string;
  pmid?: string;
  title: string;
  authors?: string[];
  journal?: string;
  year?: number;
  doi?: string;
  evidenceLevel?: string;
}

// ============================================
// Schema Definition: World 2 - Traditional Medicine
// ============================================

interface HerbNode {
  id: string;
  name: string;
  latinName?: string;
  englishName?: string;
  modality: 'tcm' | 'ayurveda' | 'kampo' | 'western' | 'other';
  tcmspId?: string;
  properties?: string[];
  meridians?: string[];
  temperature?: string;
  taste?: string;
  ayurvedaRasa?: string;
  ayurvedaVirya?: string;
  ayurvedaVipaka?: string;
}

interface FormulaNode {
  id: string;
  name: string;
  englishName?: string;
  modality: 'tcm' | 'ayurveda' | 'kampo' | 'western';
  category?: string;
  indication?: string;
  kampoId?: string;
  tcmspId?: string;
}

interface ModalityNode {
  id: string;
  name: string;
  type: 'tcm' | 'ayurveda' | 'kampo' | 'western' | 'unani' | 'siddha';
  region?: string;
  description?: string;
}

interface TCMPatternNode {
  id: string;
  name: string;
  pinyin?: string;
  category?: string;
  symptoms?: string[];
  tongue?: string;
  pulse?: string;
}

interface MeridianNode {
  id: string;
  name: string;
  pinyin?: string;
  element?: string;
  yinYang?: 'yin' | 'yang';
  pairedMeridian?: string;
}

// ============================================
// Schema Definition: World 3 - Nutrition
// ============================================

interface NutrientNode {
  id: string;
  name: string;
  alternateNames?: string[];
  category: 'vitamin' | 'mineral' | 'amino_acid' | 'fatty_acid' | 'phytonutrient' | 'other';
  subcategory?: string;
  usdaId?: string;
  recommendedDaily?: string;
  unit?: string;
}

interface FoodNode {
  id: string;
  name: string;
  category?: string;
  usdaFdcId?: string;
  foodDbId?: string;
  scientificName?: string;
}

interface BiomarkerNode {
  id: string;
  name: string;
  category: 'blood' | 'urine' | 'stool' | 'saliva' | 'genetic' | 'other';
  unit?: string;
  normalRangeLow?: number;
  normalRangeHigh?: number;
  loincCode?: string;
}

interface SupplementNode {
  id: string;
  name: string;
  type: 'vitamin' | 'mineral' | 'herb' | 'compound' | 'probiotic' | 'enzyme' | 'other';
  form?: string;
  dsldId?: string;
  standardDose?: string;
}

interface FoodGroupNode {
  id: string;
  name: string;
  description?: string;
  usdaGroupCode?: string;
}

// ============================================
// Edge Types
// ============================================

type EdgeType =
  // Genetics internal
  | 'HAS_VARIANT'        // Gene -> SNP
  | 'INVOLVED_IN'        // Gene -> Pathway
  | 'ASSOCIATED_WITH'    // SNP -> Condition
  | 'REGULATES'          // Gene -> Gene
  // Traditional internal
  | 'INGREDIENT_OF'      // Herb -> Formula
  | 'BELONGS_TO'         // Herb -> Modality
  | 'ENTERS'             // Herb -> Meridian
  | 'PATTERN_OF'         // TCMPattern -> Condition
  // Nutrition internal
  | 'CONTAINS'           // Food -> Nutrient
  | 'MEASURES'           // Biomarker -> Nutrient
  | 'IN_GROUP'           // Food -> FoodGroup
  // Cross-world
  | 'AFFECTS'            // SNP -> Nutrient (genetics -> nutrition)
  | 'TREATS'             // Herb -> Condition (traditional -> genetics)
  | 'REQUIRES_COFACTOR'  // Pathway -> Nutrient (genetics -> nutrition)
  | 'TARGETS'            // Herb -> Pathway (traditional -> genetics)
  | 'MODULATES'          // Nutrient -> Gene (nutrition -> genetics)
  | 'DEFICIENCY_CAUSES'  // Nutrient -> Condition (nutrition -> genetics)
  | 'HAS_RESEARCH'       // Any -> Publication
  | 'INTERACTS_WITH';    // Herb -> Herb, Nutrient -> Nutrient

// ============================================
// Sample Data for Schema Validation
// ============================================

const sampleGenes: GeneNode[] = [
  {
    id: 'gene_MTHFR',
    symbol: 'MTHFR',
    name: 'Methylenetetrahydrofolate reductase',
    chromosome: '1',
    description: 'Enzyme in folate metabolism and methylation cycle',
    ensemblId: 'ENSG00000177000',
    uniprotId: 'P42898',
    hgncId: 'HGNC:7436',
  },
  {
    id: 'gene_COMT',
    symbol: 'COMT',
    name: 'Catechol-O-methyltransferase',
    chromosome: '22',
    description: 'Enzyme that degrades catecholamines',
    ensemblId: 'ENSG00000093010',
    uniprotId: 'P21964',
  },
  {
    id: 'gene_VDR',
    symbol: 'VDR',
    name: 'Vitamin D receptor',
    chromosome: '12',
    description: 'Nuclear receptor for vitamin D',
    ensemblId: 'ENSG00000111424',
  },
];

const sampleSNPs: SNPNode[] = [
  {
    id: 'snp_rs1801133',
    rsid: 'rs1801133',
    geneSymbol: 'MTHFR',
    chromosome: '1',
    position: 11856378,
    refAllele: 'G',
    altAllele: 'A',
    minorAllele: 'A',
    maf: 0.34,
    clinicalSignificance: 'drug_response',
    commonName: 'MTHFR C677T',
  },
  {
    id: 'snp_rs4680',
    rsid: 'rs4680',
    geneSymbol: 'COMT',
    chromosome: '22',
    refAllele: 'G',
    altAllele: 'A',
    commonName: 'COMT Val158Met',
  },
];

const samplePathways: PathwayNode[] = [
  {
    id: 'pathway_methylation',
    name: 'Methylation Cycle',
    description: 'One-carbon metabolism and methyl group transfer',
    category: 'metabolism',
    reactomeId: 'R-HSA-196741',
    keggId: 'hsa00670',
  },
  {
    id: 'pathway_folate',
    name: 'Folate Metabolism',
    description: 'Folate biosynthesis and one-carbon pool',
    category: 'metabolism',
    keggId: 'hsa00790',
  },
];

const sampleHerbs: HerbNode[] = [
  {
    id: 'herb_huangqi',
    name: 'Huang Qi',
    latinName: 'Astragalus membranaceus',
    englishName: 'Astragalus',
    modality: 'tcm',
    tcmspId: 'MOL000239',
    properties: ['qi_tonic', 'immune_support'],
    meridians: ['spleen', 'lung'],
    temperature: 'warm',
    taste: 'sweet',
  },
  {
    id: 'herb_ashwagandha',
    name: 'Ashwagandha',
    latinName: 'Withania somnifera',
    englishName: 'Indian Ginseng',
    modality: 'ayurveda',
    properties: ['adaptogen', 'rasayana'],
    ayurvedaRasa: 'tikta_kashaya',
    ayurvedaVirya: 'ushna',
    ayurvedaVipaka: 'madhura',
  },
];

const sampleFormulas: FormulaNode[] = [
  {
    id: 'formula_buzhongyiqi',
    name: 'Bu Zhong Yi Qi Tang',
    englishName: 'Tonify the Middle and Augment the Qi Decoction',
    modality: 'tcm',
    category: 'qi_tonifying',
    indication: 'spleen_qi_deficiency',
  },
];

const sampleNutrients: NutrientNode[] = [
  {
    id: 'nutrient_methylfolate',
    name: 'Methylfolate',
    alternateNames: ['5-MTHF', 'L-methylfolate', 'Metafolin'],
    category: 'vitamin',
    subcategory: 'B_vitamin',
    usdaId: '1177',
    recommendedDaily: '400-800mcg',
    unit: 'mcg',
  },
  {
    id: 'nutrient_b12',
    name: 'Vitamin B12',
    alternateNames: ['Cobalamin', 'Methylcobalamin'],
    category: 'vitamin',
    subcategory: 'B_vitamin',
    recommendedDaily: '2.4mcg',
    unit: 'mcg',
  },
  {
    id: 'nutrient_magnesium',
    name: 'Magnesium',
    category: 'mineral',
    recommendedDaily: '310-420mg',
    unit: 'mg',
  },
];

const sampleConditions: ConditionNode[] = [
  {
    id: 'condition_mecfs',
    name: 'ME/CFS',
    fullName: 'Myalgic Encephalomyelitis / Chronic Fatigue Syndrome',
    icd10Code: 'G93.32',
    meshId: 'D015673',
    category: 'neurological',
  },
  {
    id: 'condition_homocysteinemia',
    name: 'Hyperhomocysteinemia',
    fullName: 'Elevated Homocysteine',
    icd10Code: 'E72.11',
    category: 'metabolic',
  },
];

// ============================================
// Schema Initialization
// ============================================

async function initializeSchema(): Promise<void> {
  console.log('===========================================');
  console.log('THREE WORLDS Schema Initialization');
  console.log('===========================================');
  console.log(`Storage: ${config.storagePath}`);
  console.log(`Dimensions: ${config.dimensions}`);
  console.log(`Metric: ${config.distanceMetric}`);
  console.log('');

  // Initialize database
  const db = new GraphDatabase({
    distanceMetric: config.distanceMetric,
    dimensions: config.dimensions,
    storagePath: config.storagePath,
  });

  console.log(`Database initialized. Persistent: ${db.isPersistent()}`);
  console.log('');

  // ==========================================
  // World 1: Genetics
  // ==========================================
  console.log('--- World 1: Genetics ---');

  // Create Gene nodes
  for (const gene of sampleGenes) {
    await db.createNode({
      id: gene.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['Gene'],
      properties: {
        symbol: gene.symbol,
        name: gene.name,
        chromosome: gene.chromosome || '',
        description: gene.description || '',
        ensemblId: gene.ensemblId || '',
        uniprotId: gene.uniprotId || '',
        hgncId: gene.hgncId || '',
      },
    });
    console.log(`  Created Gene: ${gene.symbol}`);
  }

  // Create SNP nodes
  for (const snp of sampleSNPs) {
    await db.createNode({
      id: snp.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['SNP'],
      properties: {
        rsid: snp.rsid,
        geneSymbol: snp.geneSymbol || '',
        chromosome: snp.chromosome || '',
        position: String(snp.position || 0),
        refAllele: snp.refAllele || '',
        altAllele: snp.altAllele || '',
        clinicalSignificance: snp.clinicalSignificance || '',
        commonName: snp.commonName || '',
      },
    });
    console.log(`  Created SNP: ${snp.rsid}`);
  }

  // Create Pathway nodes
  for (const pathway of samplePathways) {
    await db.createNode({
      id: pathway.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['Pathway'],
      properties: {
        name: pathway.name,
        description: pathway.description || '',
        category: pathway.category || '',
        reactomeId: pathway.reactomeId || '',
        keggId: pathway.keggId || '',
      },
    });
    console.log(`  Created Pathway: ${pathway.name}`);
  }

  // Create Condition nodes
  for (const condition of sampleConditions) {
    await db.createNode({
      id: condition.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['Condition'],
      properties: {
        name: condition.name,
        fullName: condition.fullName || '',
        icd10Code: condition.icd10Code || '',
        meshId: condition.meshId || '',
        category: condition.category || '',
      },
    });
    console.log(`  Created Condition: ${condition.name}`);
  }

  // ==========================================
  // World 2: Traditional Medicine
  // ==========================================
  console.log('');
  console.log('--- World 2: Traditional Medicine ---');

  // Create Herb nodes
  for (const herb of sampleHerbs) {
    await db.createNode({
      id: herb.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['Herb'],
      properties: {
        name: herb.name,
        latinName: herb.latinName || '',
        englishName: herb.englishName || '',
        modality: herb.modality,
        tcmspId: herb.tcmspId || '',
        properties: JSON.stringify(herb.properties || []),
        meridians: JSON.stringify(herb.meridians || []),
        temperature: herb.temperature || '',
        taste: herb.taste || '',
      },
    });
    console.log(`  Created Herb: ${herb.name} (${herb.modality})`);
  }

  // Create Formula nodes
  for (const formula of sampleFormulas) {
    await db.createNode({
      id: formula.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['Formula'],
      properties: {
        name: formula.name,
        englishName: formula.englishName || '',
        modality: formula.modality,
        category: formula.category || '',
        indication: formula.indication || '',
      },
    });
    console.log(`  Created Formula: ${formula.name}`);
  }

  // ==========================================
  // World 3: Nutrition
  // ==========================================
  console.log('');
  console.log('--- World 3: Nutrition ---');

  // Create Nutrient nodes
  for (const nutrient of sampleNutrients) {
    await db.createNode({
      id: nutrient.id,
      embedding: randomEmbedding(config.dimensions),
      labels: ['Nutrient'],
      properties: {
        name: nutrient.name,
        alternateNames: JSON.stringify(nutrient.alternateNames || []),
        category: nutrient.category,
        subcategory: nutrient.subcategory || '',
        usdaId: nutrient.usdaId || '',
        recommendedDaily: nutrient.recommendedDaily || '',
        unit: nutrient.unit || '',
      },
    });
    console.log(`  Created Nutrient: ${nutrient.name}`);
  }

  // ==========================================
  // Create Edges (Relationships)
  // ==========================================
  console.log('');
  console.log('--- Creating Relationships ---');

  // Gene -> SNP (HAS_VARIANT)
  await db.createEdge({
    from: 'gene_MTHFR',
    to: 'snp_rs1801133',
    description: 'HAS_VARIANT',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
    metadata: { source: 'dbSNP' },
  });
  console.log('  Created: MTHFR -[HAS_VARIANT]-> rs1801133');

  await db.createEdge({
    from: 'gene_COMT',
    to: 'snp_rs4680',
    description: 'HAS_VARIANT',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
  });
  console.log('  Created: COMT -[HAS_VARIANT]-> rs4680');

  // Gene -> Pathway (INVOLVED_IN)
  await db.createEdge({
    from: 'gene_MTHFR',
    to: 'pathway_methylation',
    description: 'INVOLVED_IN',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
    metadata: { role: 'enzyme' },
  });
  console.log('  Created: MTHFR -[INVOLVED_IN]-> Methylation Cycle');

  await db.createEdge({
    from: 'gene_MTHFR',
    to: 'pathway_folate',
    description: 'INVOLVED_IN',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
  });
  console.log('  Created: MTHFR -[INVOLVED_IN]-> Folate Metabolism');

  // SNP -> Nutrient (AFFECTS) - Cross-world!
  await db.createEdge({
    from: 'snp_rs1801133',
    to: 'nutrient_methylfolate',
    description: 'AFFECTS',
    embedding: randomEmbedding(config.dimensions),
    confidence: 0.95,
    metadata: { effect: 'reduces_conversion', magnitude: 'moderate' },
  });
  console.log('  Created: rs1801133 -[AFFECTS]-> Methylfolate (Cross-world: Genetics->Nutrition)');

  // Pathway -> Nutrient (REQUIRES_COFACTOR) - Cross-world!
  await db.createEdge({
    from: 'pathway_methylation',
    to: 'nutrient_b12',
    description: 'REQUIRES_COFACTOR',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
  });
  console.log('  Created: Methylation -[REQUIRES_COFACTOR]-> B12 (Cross-world: Genetics->Nutrition)');

  await db.createEdge({
    from: 'pathway_methylation',
    to: 'nutrient_methylfolate',
    description: 'REQUIRES_COFACTOR',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
  });
  console.log('  Created: Methylation -[REQUIRES_COFACTOR]-> Methylfolate');

  // SNP -> Condition (ASSOCIATED_WITH)
  await db.createEdge({
    from: 'snp_rs1801133',
    to: 'condition_homocysteinemia',
    description: 'ASSOCIATED_WITH',
    embedding: randomEmbedding(config.dimensions),
    confidence: 0.85,
    metadata: { oddsRatio: '1.5', pValue: '0.001' },
  });
  console.log('  Created: rs1801133 -[ASSOCIATED_WITH]-> Hyperhomocysteinemia');

  // Herb -> Formula (INGREDIENT_OF)
  await db.createEdge({
    from: 'herb_huangqi',
    to: 'formula_buzhongyiqi',
    description: 'INGREDIENT_OF',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
    metadata: { role: 'chief', amount: '15g' },
  });
  console.log('  Created: Huang Qi -[INGREDIENT_OF]-> Bu Zhong Yi Qi Tang');

  // Herb -> Condition (TREATS) - Cross-world!
  await db.createEdge({
    from: 'herb_ashwagandha',
    to: 'condition_mecfs',
    description: 'TREATS',
    embedding: randomEmbedding(config.dimensions),
    confidence: 0.7,
    metadata: { evidenceLevel: 'clinical', mechanism: 'adaptogenic' },
  });
  console.log('  Created: Ashwagandha -[TREATS]-> ME/CFS (Cross-world: Traditional->Genetics)');

  // ==========================================
  // Create Hyperedges (Multi-node relationships)
  // ==========================================
  console.log('');
  console.log('--- Creating Hyperedges ---');

  // Formula composition hyperedge
  await db.createHyperedge({
    nodes: ['herb_huangqi', 'formula_buzhongyiqi'],
    description: 'FORMULA_COMPOSITION',
    embedding: randomEmbedding(config.dimensions),
    confidence: 1.0,
    metadata: { formulaType: 'qi_tonifying' },
  });
  console.log('  Created Hyperedge: Formula Composition (Huang Qi + Bu Zhong Yi Qi Tang)');

  // Gene-Nutrient-Condition triangle
  await db.createHyperedge({
    nodes: ['gene_MTHFR', 'nutrient_methylfolate', 'condition_homocysteinemia'],
    description: 'GENE_NUTRIENT_CONDITION_PATHWAY',
    embedding: randomEmbedding(config.dimensions),
    confidence: 0.9,
    metadata: {
      mechanism: 'MTHFR variants reduce methylfolate conversion, leading to elevated homocysteine',
    },
  });
  console.log('  Created Hyperedge: Gene-Nutrient-Condition (MTHFR + Methylfolate + Homocysteinemia)');

  // ==========================================
  // Verify Schema
  // ==========================================
  console.log('');
  console.log('--- Verification ---');

  const stats = await db.stats();
  console.log(`Total Nodes: ${stats.totalNodes}`);
  console.log(`Total Edges: ${stats.totalEdges}`);
  console.log(`Average Degree: ${stats.avgDegree.toFixed(2)}`);

  // Test Cypher query
  console.log('');
  console.log('--- Test Queries ---');

  const geneResult = await db.query('MATCH (g:Gene) RETURN g');
  console.log(`Genes found: ${geneResult.nodes.length}`);

  const snpResult = await db.query('MATCH (s:SNP) RETURN s');
  console.log(`SNPs found: ${snpResult.nodes.length}`);

  const herbResult = await db.query('MATCH (h:Herb) RETURN h');
  console.log(`Herbs found: ${herbResult.nodes.length}`);

  const nutrientResult = await db.query('MATCH (n:Nutrient) RETURN n');
  console.log(`Nutrients found: ${nutrientResult.nodes.length}`);

  // Test cross-world query
  const crossWorldResult = await db.query(`
    MATCH (s:SNP)-[:AFFECTS]->(n:Nutrient)
    RETURN s, n
  `);
  console.log(`Cross-world SNP-Nutrient relationships: ${crossWorldResult.edges.length}`);

  console.log('');
  console.log('===========================================');
  console.log('Schema initialization complete!');
  console.log('===========================================');
}

// Run initialization
initializeSchema()
  .then(() => process.exit(0))
  .catch((error) => {
    console.error('Schema initialization failed:', error);
    process.exit(1);
  });
