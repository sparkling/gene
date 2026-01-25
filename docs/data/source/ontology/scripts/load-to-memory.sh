#!/bin/bash
# Load knowledge graph summaries into claude-flow memory for development use
# Run from project root: bash docs/data/source/ontology/scripts/load-to-memory.sh

set -e
cd "$(dirname "$0")/../../../.."

echo "=== Loading Knowledge Graph Summaries into Claude-Flow Memory ==="

# Initialize memory if needed
npx @claude-flow/cli@latest memory init --force 2>/dev/null || true

# Store category summaries
echo "Loading category summaries..."

npx @claude-flow/cli@latest memory store \
  --key "kg-categories" \
  --value "9 categories: genetics-genomics, compounds-molecules, diseases-phenotypes, pathways-networks, traditional-medicine, nutrition-food, proteins-molecular, literature-knowledge, microbiome" \
  --namespace knowledge

npx @claude-flow/cli@latest memory store \
  --key "kg-statistics" \
  --value "137 DataSources, 103 Licenses, 58 Versions, 257 AccessMethods, 232 Identifiers, 134 CrossReferences, 324 UseCases, 129 Limitations across 10,832 lines of Turtle RDF" \
  --namespace knowledge

# Store tier information
echo "Loading tier information..."

npx @claude-flow/cli@latest memory store \
  --key "tier1-sources" \
  --value "Tier 1 (Essential): ClinVar, gnomAD, dbSNP, Ensembl, UniProt, KEGG, Reactome, PubChem, DrugBank, ChEMBL, OMIM, HPO, DisGeNET, PubMed, Gene Ontology" \
  --namespace knowledge

npx @claude-flow/cli@latest memory store \
  --key "tier2-sources" \
  --value "Tier 2 (Important): CADD, REVEL, GTEx, PharmGKB, TCMSP, SymMap, FooDB, USDA, STRING, BioGRID" \
  --namespace knowledge

# Store size information
echo "Loading data size information..."

npx @claude-flow/cli@latest memory store \
  --key "large-datasets" \
  --value "gnomAD: 1.5TB, dbSNP: 50GB, CADD: 85GB, UniProt: 5GB, PubChem: 50GB, ChEMBL: 10GB, Ensembl: 40GB zipped" \
  --namespace knowledge

npx @claude-flow/cli@latest memory store \
  --key "medium-datasets" \
  --value "ClinVar: 2.5GB, KEGG: 1GB, Reactome: 500MB, DrugBank: 150MB, OMIM: 200MB, HPO: 100MB" \
  --namespace knowledge

# Store access methods
echo "Loading access method information..."

npx @claude-flow/cli@latest memory store \
  --key "access-rest-api" \
  --value "REST API sources: ClinVar (E-utils), gnomAD (GraphQL), Ensembl (REST), UniProt (REST), KEGG (REST), PubChem (PUG-REST), DrugBank (REST), OMIM (REST), PubMed (E-utils)" \
  --namespace knowledge

npx @claude-flow/cli@latest memory store \
  --key "access-ftp" \
  --value "FTP download: ClinVar, gnomAD, dbSNP, Ensembl, UniProt, KEGG, Reactome, PubChem, ChEMBL, GTEx" \
  --namespace knowledge

npx @claude-flow/cli@latest memory store \
  --key "access-sparql" \
  --value "SPARQL endpoints: UniProt, Wikidata, DisGeNET, Gene Ontology, ChEBI" \
  --namespace knowledge

# Store cross-reference information
echo "Loading cross-reference patterns..."

npx @claude-flow/cli@latest memory store \
  --key "xref-patterns" \
  --value "Common cross-references: rsID (dbSNP), HGVS (variants), Ensembl ID (genes), UniProt ID (proteins), KEGG ID (pathways), PubChem CID (compounds), ChEMBL ID (drugs), OMIM ID (diseases), HPO ID (phenotypes)" \
  --namespace knowledge

# Store ontology file locations
echo "Loading file locations..."

npx @claude-flow/cli@latest memory store \
  --key "ontology-files" \
  --value "OWL: docs/data/source/ontology/owl/datasource-ontology.ttl | SHACL: docs/data/source/ontology/shacl/datasource-shapes.ttl | SKOS: docs/data/source/ontology/skos/datasource-taxonomy.ttl | Instances: docs/data/source/ontology/instances/0*.ttl" \
  --namespace knowledge

# Store SPARQL query examples
echo "Loading query examples..."

npx @claude-flow/cli@latest memory store \
  --key "sparql-examples" \
  --value "Find sources with sizes: SELECT ?source ?size WHERE { ?v a ds:Version ; ds:forSource ?source ; ds:totalSize ?size } | Find by tier: SELECT ?source WHERE { ?source ds:tier 1 } | Find by access: SELECT ?source ?url WHERE { ?a ds:forSource ?source ; ds:methodType 'REST API' ; ds:baseUrl ?url }" \
  --namespace knowledge

echo ""
echo "=== Knowledge Graph loaded into claude-flow memory ==="
echo ""
echo "Query examples:"
echo "  npx @claude-flow/cli@latest memory search --query 'genetic variant databases'"
echo "  npx @claude-flow/cli@latest memory search --query 'large datasets over 1TB'"
echo "  npx @claude-flow/cli@latest memory search --query 'REST API access methods'"
echo "  npx @claude-flow/cli@latest memory retrieve --key 'tier1-sources' --namespace knowledge"
