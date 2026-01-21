/**
 * RuVector Graph Demo for Gene Project
 *
 * Demonstrates the complete workflow:
 * 1. Initialize graph database
 * 2. Create nodes (Gene, SNP, Pathway, Nutrient, Herb)
 * 3. Create edges with properties
 * 4. Create hyperedges for complex relationships
 * 5. Run Cypher queries
 * 6. Perform vector similarity search
 * 7. Transaction management
 *
 * Usage:
 *   npx ts-node src/lib/ruvector-graph-demo.ts
 *
 * @module ruvector-graph-demo
 */

import { RuVectorGraph, randomEmbedding } from './ruvector-graph';
import { RuVectorEmbeddings } from './ruvector-embeddings';

// ============================================
// Demo Configuration
// ============================================

const DEMO_CONFIG = {
  storagePath: './data/demo-knowledge.db',
  dimensions: 384,
  distanceMetric: 'Cosine' as const,
};

// ============================================
// Demo Data: THREE WORLDS
// ============================================

/**
 * World 1: Genetics sample data
 */
const geneticsData = {
  genes: [
    {
      id: 'gene_MTHFR',
      symbol: 'MTHFR',
      name: 'Methylenetetrahydrofolate reductase',
      chromosome: '1',
      description: 'Key enzyme in folate metabolism and methylation',
    },
    {
      id: 'gene_COMT',
      symbol: 'COMT',
      name: 'Catechol-O-methyltransferase',
      chromosome: '22',
      description: 'Degrades catecholamines including dopamine',
    },
    {
      id: 'gene_VDR',
      symbol: 'VDR',
      name: 'Vitamin D receptor',
      chromosome: '12',
      description: 'Nuclear receptor that mediates vitamin D effects',
    },
  ],
  snps: [
    {
      id: 'snp_rs1801133',
      rsid: 'rs1801133',
      geneSymbol: 'MTHFR',
      commonName: 'MTHFR C677T',
      clinicalSignificance: 'drug_response',
    },
    {
      id: 'snp_rs4680',
      rsid: 'rs4680',
      geneSymbol: 'COMT',
      commonName: 'COMT Val158Met',
      clinicalSignificance: 'drug_response',
    },
  ],
  pathways: [
    {
      id: 'pathway_methylation',
      name: 'Methylation Cycle',
      description: 'One-carbon metabolism and methyl transfer',
      category: 'metabolism',
    },
    {
      id: 'pathway_folate',
      name: 'Folate Metabolism',
      description: 'Folate biosynthesis and utilization',
      category: 'metabolism',
    },
  ],
};

/**
 * World 2: Traditional Medicine sample data
 */
const traditionalData = {
  herbs: [
    {
      id: 'herb_huangqi',
      name: 'Huang Qi',
      latinName: 'Astragalus membranaceus',
      englishName: 'Astragalus',
      modality: 'tcm',
      properties: ['qi_tonic', 'immune_support'],
    },
    {
      id: 'herb_ashwagandha',
      name: 'Ashwagandha',
      latinName: 'Withania somnifera',
      modality: 'ayurveda',
      properties: ['adaptogen', 'rasayana'],
    },
  ],
  formulas: [
    {
      id: 'formula_buzhongyiqi',
      name: 'Bu Zhong Yi Qi Tang',
      englishName: 'Tonify Middle Qi Decoction',
      modality: 'tcm',
      indication: 'spleen_qi_deficiency',
    },
  ],
};

/**
 * World 3: Nutrition sample data
 */
const nutritionData = {
  nutrients: [
    {
      id: 'nutrient_methylfolate',
      name: 'Methylfolate',
      alternateNames: ['5-MTHF', 'L-methylfolate'],
      category: 'vitamin',
      subcategory: 'B_vitamin',
    },
    {
      id: 'nutrient_b12',
      name: 'Vitamin B12',
      alternateNames: ['Cobalamin', 'Methylcobalamin'],
      category: 'vitamin',
      subcategory: 'B_vitamin',
    },
    {
      id: 'nutrient_magnesium',
      name: 'Magnesium',
      category: 'mineral',
    },
  ],
  conditions: [
    {
      id: 'condition_homocysteinemia',
      name: 'Hyperhomocysteinemia',
      fullName: 'Elevated Homocysteine',
      category: 'metabolic',
    },
    {
      id: 'condition_mecfs',
      name: 'ME/CFS',
      fullName: 'Myalgic Encephalomyelitis / Chronic Fatigue Syndrome',
      category: 'neurological',
    },
  ],
};

// ============================================
// Demo Functions
// ============================================

/**
 * Create nodes from demo data
 */
async function createNodes(
  graph: RuVectorGraph,
  embeddings: RuVectorEmbeddings
): Promise<void> {
  console.log('\n--- Creating Nodes ---');

  // World 1: Genetics
  console.log('\nWorld 1: Genetics');

  for (const gene of geneticsData.genes) {
    const embedding = await embeddings.embedGene(gene);
    await graph.createNode({
      id: gene.id,
      embedding,
      labels: ['Gene'],
      properties: {
        symbol: gene.symbol,
        name: gene.name,
        chromosome: gene.chromosome,
        description: gene.description,
      },
    });
    console.log(`  Created Gene: ${gene.symbol}`);
  }

  for (const snp of geneticsData.snps) {
    const embedding = await embeddings.embedSNP(snp);
    await graph.createNode({
      id: snp.id,
      embedding,
      labels: ['SNP'],
      properties: {
        rsid: snp.rsid,
        geneSymbol: snp.geneSymbol,
        commonName: snp.commonName,
        clinicalSignificance: snp.clinicalSignificance,
      },
    });
    console.log(`  Created SNP: ${snp.rsid}`);
  }

  for (const pathway of geneticsData.pathways) {
    const embedding = await embeddings.embedPathway(pathway);
    await graph.createNode({
      id: pathway.id,
      embedding,
      labels: ['Pathway'],
      properties: {
        name: pathway.name,
        description: pathway.description,
        category: pathway.category,
      },
    });
    console.log(`  Created Pathway: ${pathway.name}`);
  }

  // World 2: Traditional Medicine
  console.log('\nWorld 2: Traditional Medicine');

  for (const herb of traditionalData.herbs) {
    const embedding = await embeddings.embedHerb(herb);
    await graph.createNode({
      id: herb.id,
      embedding,
      labels: ['Herb'],
      properties: {
        name: herb.name,
        latinName: herb.latinName || '',
        englishName: herb.englishName || '',
        modality: herb.modality,
        properties: JSON.stringify(herb.properties),
      },
    });
    console.log(`  Created Herb: ${herb.name} (${herb.modality})`);
  }

  for (const formula of traditionalData.formulas) {
    const embedding = await embeddings.embed(
      `Formula: ${formula.name}. ${formula.englishName}. Indication: ${formula.indication}`
    );
    await graph.createNode({
      id: formula.id,
      embedding,
      labels: ['Formula'],
      properties: {
        name: formula.name,
        englishName: formula.englishName,
        modality: formula.modality,
        indication: formula.indication,
      },
    });
    console.log(`  Created Formula: ${formula.name}`);
  }

  // World 3: Nutrition
  console.log('\nWorld 3: Nutrition');

  for (const nutrient of nutritionData.nutrients) {
    const embedding = await embeddings.embedNutrient(nutrient);
    await graph.createNode({
      id: nutrient.id,
      embedding,
      labels: ['Nutrient'],
      properties: {
        name: nutrient.name,
        alternateNames: JSON.stringify(nutrient.alternateNames || []),
        category: nutrient.category,
        subcategory: nutrient.subcategory || '',
      },
    });
    console.log(`  Created Nutrient: ${nutrient.name}`);
  }

  for (const condition of nutritionData.conditions) {
    const embedding = await embeddings.embed(
      `Condition: ${condition.name}. ${condition.fullName}. Category: ${condition.category}`
    );
    await graph.createNode({
      id: condition.id,
      embedding,
      labels: ['Condition'],
      properties: {
        name: condition.name,
        fullName: condition.fullName,
        category: condition.category,
      },
    });
    console.log(`  Created Condition: ${condition.name}`);
  }
}

/**
 * Create edges between nodes
 */
async function createEdges(
  graph: RuVectorGraph,
  embeddings: RuVectorEmbeddings
): Promise<void> {
  console.log('\n--- Creating Edges ---');

  const edges = [
    // Genetics internal: Gene -> SNP
    {
      from: 'gene_MTHFR',
      to: 'snp_rs1801133',
      type: 'HAS_VARIANT',
      metadata: { source: 'dbSNP' },
    },
    {
      from: 'gene_COMT',
      to: 'snp_rs4680',
      type: 'HAS_VARIANT',
      metadata: { source: 'dbSNP' },
    },
    // Genetics internal: Gene -> Pathway
    {
      from: 'gene_MTHFR',
      to: 'pathway_methylation',
      type: 'INVOLVED_IN',
      metadata: { role: 'enzyme' },
    },
    {
      from: 'gene_MTHFR',
      to: 'pathway_folate',
      type: 'INVOLVED_IN',
      metadata: { role: 'catalyst' },
    },
    // Cross-world: SNP -> Nutrient (Genetics -> Nutrition)
    {
      from: 'snp_rs1801133',
      to: 'nutrient_methylfolate',
      type: 'AFFECTS',
      metadata: { effect: 'reduces_conversion', magnitude: 'moderate' },
      confidence: 0.95,
    },
    // Cross-world: Pathway -> Nutrient (Genetics -> Nutrition)
    {
      from: 'pathway_methylation',
      to: 'nutrient_b12',
      type: 'REQUIRES_COFACTOR',
      metadata: {},
    },
    {
      from: 'pathway_methylation',
      to: 'nutrient_methylfolate',
      type: 'REQUIRES_COFACTOR',
      metadata: {},
    },
    // Genetics internal: SNP -> Condition
    {
      from: 'snp_rs1801133',
      to: 'condition_homocysteinemia',
      type: 'ASSOCIATED_WITH',
      metadata: { oddsRatio: 1.5, pValue: 0.001 },
      confidence: 0.85,
    },
    // Traditional internal: Herb -> Formula
    {
      from: 'herb_huangqi',
      to: 'formula_buzhongyiqi',
      type: 'INGREDIENT_OF',
      metadata: { role: 'chief', amount: '15g' },
    },
    // Cross-world: Herb -> Condition (Traditional -> Genetics)
    {
      from: 'herb_ashwagandha',
      to: 'condition_mecfs',
      type: 'TREATS',
      metadata: { evidenceLevel: 'clinical', mechanism: 'adaptogenic' },
      confidence: 0.7,
    },
  ];

  for (const edge of edges) {
    const relEmbed = await embeddings.embedRelationship({
      type: edge.type,
      fromLabel: edge.from.split('_')[0],
      toLabel: edge.to.split('_')[0],
      metadata: edge.metadata,
    });

    await graph.createEdge({
      from: edge.from,
      to: edge.to,
      description: edge.type,
      embedding: relEmbed,
      confidence: edge.confidence || 1.0,
      metadata: edge.metadata,
    });

    console.log(`  Created: ${edge.from} -[${edge.type}]-> ${edge.to}`);
  }
}

/**
 * Create hyperedges for complex relationships
 */
async function createHyperedges(
  graph: RuVectorGraph,
  embeddings: RuVectorEmbeddings
): Promise<void> {
  console.log('\n--- Creating Hyperedges ---');

  // Gene-Nutrient-Condition triangle
  const geneNutrientCondition = await embeddings.embed(
    'Gene nutrient condition pathway: MTHFR variants reduce methylfolate conversion, leading to elevated homocysteine'
  );

  await graph.createHyperedge({
    nodes: ['gene_MTHFR', 'nutrient_methylfolate', 'condition_homocysteinemia'],
    description: 'GENE_NUTRIENT_CONDITION_PATHWAY',
    embedding: geneNutrientCondition,
    confidence: 0.9,
    metadata: {
      mechanism: 'MTHFR variants reduce methylfolate conversion',
      clinicalRelevance: 'high',
    },
  });
  console.log(
    '  Created: Gene-Nutrient-Condition hyperedge (MTHFR + Methylfolate + Homocysteinemia)'
  );

  // Formula composition hyperedge
  const formulaComposition = await embeddings.embed(
    'Traditional Chinese Medicine formula composition for tonifying middle qi'
  );

  await graph.createHyperedge({
    nodes: ['herb_huangqi', 'formula_buzhongyiqi'],
    description: 'FORMULA_COMPOSITION',
    embedding: formulaComposition,
    confidence: 1.0,
    metadata: {
      formulaType: 'qi_tonifying',
      tradition: 'tcm',
    },
  });
  console.log('  Created: Formula composition hyperedge (Huang Qi + Bu Zhong Yi Qi Tang)');

  // Cross-world nutritional genomics hyperedge
  const nutrigenomics = await embeddings.embed(
    'Nutritional genomics: genetic variants affecting nutrient metabolism and requirements'
  );

  await graph.createHyperedge({
    nodes: ['gene_MTHFR', 'snp_rs1801133', 'nutrient_b12', 'nutrient_methylfolate'],
    description: 'NUTRIGENOMICS_PATHWAY',
    embedding: nutrigenomics,
    confidence: 0.88,
    metadata: {
      field: 'nutrigenomics',
      recommendation: 'Consider methylated B vitamins for MTHFR C677T carriers',
    },
  });
  console.log(
    '  Created: Nutrigenomics hyperedge (MTHFR + rs1801133 + B12 + Methylfolate)'
  );
}

/**
 * Run Cypher queries
 */
async function runQueries(graph: RuVectorGraph): Promise<void> {
  console.log('\n--- Running Cypher Queries ---');

  // Query 1: Find all genes
  console.log('\n1. Find all genes:');
  const genes = await graph.query('MATCH (g:Gene) RETURN g');
  console.log(`   Found ${genes.nodes.length} genes`);
  genes.nodes.forEach((node) => {
    console.log(`   - ${node.properties.symbol}: ${node.properties.name}`);
  });

  // Query 2: Find all SNPs
  console.log('\n2. Find all SNPs:');
  const snps = await graph.query('MATCH (s:SNP) RETURN s');
  console.log(`   Found ${snps.nodes.length} SNPs`);
  snps.nodes.forEach((node) => {
    console.log(`   - ${node.properties.rsid} (${node.properties.commonName})`);
  });

  // Query 3: Find SNPs affecting nutrients (cross-world)
  console.log('\n3. Find SNPs affecting nutrients (cross-world query):');
  const snpNutrients = await graph.query(`
    MATCH (s:SNP)-[a:AFFECTS]->(n:Nutrient)
    RETURN s, a, n
  `);
  console.log(`   Found ${snpNutrients.edges.length} SNP-Nutrient relationships`);

  // Query 4: Find herbs and their formulas
  console.log('\n4. Find herbs in formulas:');
  const herbFormulas = await graph.query(`
    MATCH (h:Herb)-[i:INGREDIENT_OF]->(f:Formula)
    RETURN h, i, f
  `);
  console.log(`   Found ${herbFormulas.edges.length} herb-formula relationships`);

  // Query 5: Find all cross-world relationships
  console.log('\n5. Find pathway cofactors:');
  const cofactors = await graph.query(`
    MATCH (p:Pathway)-[r:REQUIRES_COFACTOR]->(n:Nutrient)
    RETURN p, r, n
  `);
  console.log(`   Found ${cofactors.edges.length} pathway-nutrient cofactor relationships`);
}

/**
 * Demonstrate graph traversal
 */
async function demonstrateTraversal(graph: RuVectorGraph): Promise<void> {
  console.log('\n--- Graph Traversal ---');

  // 1-hop neighbors of MTHFR
  console.log('\n1. Direct neighbors of MTHFR gene:');
  const neighbors1 = await graph.getNeighbors('gene_MTHFR', 1);
  console.log(`   Found ${neighbors1.length} direct neighbors`);
  neighbors1.forEach((n) => console.log(`   - ${n}`));

  // 2-hop neighbors of MTHFR
  console.log('\n2. 2-hop neighborhood of MTHFR gene:');
  const neighbors2 = await graph.getNeighbors('gene_MTHFR', 2);
  console.log(`   Found ${neighbors2.length} nodes within 2 hops`);
}

/**
 * Demonstrate vector similarity search
 */
async function demonstrateSimilaritySearch(
  graph: RuVectorGraph,
  embeddings: RuVectorEmbeddings
): Promise<void> {
  console.log('\n--- Vector Similarity Search ---');

  // Search for hyperedges similar to "methylation defects"
  console.log('\n1. Searching for hyperedges related to "methylation defects":');
  const queryEmbed = await embeddings.embed(
    'methylation defects and folate metabolism problems'
  );
  const similar = await graph.searchHyperedges(queryEmbed, 5);
  console.log(`   Found ${similar.length} similar hyperedges`);
  similar.forEach((result) => {
    console.log(`   - ${result.description || result.id} (score: ${result.score.toFixed(3)})`);
  });

  // Demonstrate embedding similarity
  console.log('\n2. Comparing embedding similarity:');
  const mthfrEmbed = await embeddings.embedGene({
    symbol: 'MTHFR',
    name: 'Methylenetetrahydrofolate reductase',
    description: 'Key enzyme in folate metabolism',
  });

  const folateEmbed = await embeddings.embedNutrient({
    name: 'Methylfolate',
    category: 'vitamin',
    alternateNames: ['5-MTHF'],
  });

  const magnesiumEmbed = await embeddings.embedNutrient({
    name: 'Magnesium',
    category: 'mineral',
  });

  const sim1 = embeddings.similarity(mthfrEmbed, folateEmbed);
  const sim2 = embeddings.similarity(mthfrEmbed, magnesiumEmbed);

  console.log(`   MTHFR <-> Methylfolate similarity: ${sim1.toFixed(3)}`);
  console.log(`   MTHFR <-> Magnesium similarity: ${sim2.toFixed(3)}`);
  console.log(`   (Higher value = more semantically related)`);
}

/**
 * Demonstrate transaction management
 */
async function demonstrateTransactions(graph: RuVectorGraph): Promise<void> {
  console.log('\n--- Transaction Management ---');

  // Successful transaction
  console.log('\n1. Successful transaction:');
  try {
    await graph.withTransaction(async (txId) => {
      console.log(`   Started transaction: ${txId}`);
      // Operations would go here
      console.log('   Committing transaction...');
    });
    console.log('   Transaction committed successfully');
  } catch (error) {
    console.log(`   Transaction failed: ${error}`);
  }

  // Rollback demonstration (simulated)
  console.log('\n2. Rollback scenario:');
  try {
    await graph.withTransaction(async (txId) => {
      console.log(`   Started transaction: ${txId}`);
      // Simulate an error
      throw new Error('Simulated error to trigger rollback');
    });
  } catch (error) {
    console.log(`   Transaction rolled back due to: ${(error as Error).message}`);
  }
}

/**
 * Display graph statistics
 */
async function displayStats(graph: RuVectorGraph): Promise<void> {
  console.log('\n--- Graph Statistics ---');

  const stats = await graph.getStats();
  console.log(`   Total Nodes: ${stats.totalNodes}`);
  console.log(`   Total Edges: ${stats.totalEdges}`);
  console.log(`   Total Hyperedges: ${stats.totalHyperedges}`);
  console.log(`   Average Degree: ${stats.avgDegree.toFixed(2)}`);
  console.log(`   Persistent: ${graph.isPersistent()}`);

  const storagePath = graph.getStoragePath();
  if (storagePath) {
    console.log(`   Storage Path: ${storagePath}`);
  }
}

// ============================================
// Main Demo Entry Point
// ============================================

async function main(): Promise<void> {
  console.log('================================================');
  console.log('RuVector Graph Demo - THREE WORLDS Knowledge Graph');
  console.log('================================================');
  console.log(`Storage: ${DEMO_CONFIG.storagePath}`);
  console.log(`Dimensions: ${DEMO_CONFIG.dimensions}`);
  console.log(`Distance Metric: ${DEMO_CONFIG.distanceMetric}`);

  try {
    // Initialize graph database
    console.log('\n--- Initializing ---');
    const graph = new RuVectorGraph(DEMO_CONFIG);
    await graph.initialize();
    console.log('Graph database initialized');

    // Initialize embeddings
    const embeddings = new RuVectorEmbeddings();
    await embeddings.initialize();
    console.log('Embeddings system initialized');

    // Run demo steps
    await createNodes(graph, embeddings);
    await createEdges(graph, embeddings);
    await createHyperedges(graph, embeddings);
    await runQueries(graph);
    await demonstrateTraversal(graph);
    await demonstrateSimilaritySearch(graph, embeddings);
    await demonstrateTransactions(graph);
    await displayStats(graph);

    console.log('\n================================================');
    console.log('Demo completed successfully!');
    console.log('================================================');
  } catch (error) {
    console.error('\nDemo failed:', error);
    process.exit(1);
  }
}

// Run if executed directly
if (require.main === module) {
  main().catch(console.error);
}

export { main as runDemo };
