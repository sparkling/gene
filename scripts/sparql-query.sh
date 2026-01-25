#!/bin/bash
# SPARQL Query Helper for Gene Knowledge Graph
# Usage: ./scripts/sparql-query.sh "SELECT * WHERE { ?s ?p ?o } LIMIT 10"
#
# Server runs on localhost:3030 (not exposed to internet)
# Start server: ./scripts/start-fuseki.sh

ENDPOINT="http://localhost:3030/gene/sparql"

if [ -z "$1" ]; then
    echo "Usage: $0 'SPARQL QUERY'"
    echo ""
    echo "Example queries:"
    echo "  $0 'SELECT (COUNT(*) AS ?c) WHERE { ?s ?p ?o }'"
    echo "  $0 'PREFIX ds: <https://gene.ai/ontology/datasource#>"
    echo "       SELECT ?source ?label WHERE { ?source a ds:DataSource ; rdfs:label ?label } LIMIT 10'"
    exit 1
fi

QUERY="$1"
FORMAT="${2:-json}"

curl -s "$ENDPOINT" \
    -H "Accept: application/${FORMAT}" \
    --data-urlencode "query=${QUERY}"
