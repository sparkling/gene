#!/bin/bash
# Start Apache Jena Fuseki with Gene Knowledge Graph
# Server binds to localhost only (not exposed to internet)
#
# Endpoint: http://localhost:3030/gene/sparql
# UI (optional): http://localhost:3030

FUSEKI_HOME="${HOME}/.local/jena/apache-jena-fuseki-5.3.0"
DATA_FILE="$(dirname "$0")/../data/qlever/gene-knowledge-graph.ttl"
LOG_FILE="/tmp/fuseki-gene.log"
PID_FILE="/tmp/fuseki-gene.pid"

# Check if already running
if [ -f "$PID_FILE" ]; then
    PID=$(cat "$PID_FILE")
    if kill -0 "$PID" 2>/dev/null; then
        echo "Fuseki already running (PID: $PID)"
        echo "Endpoint: http://localhost:3030/gene/sparql"
        exit 0
    fi
fi

# Check data file exists
if [ ! -f "$DATA_FILE" ]; then
    echo "Error: Data file not found: $DATA_FILE"
    echo "Run: source ~/.venv/qlever/bin/activate && python3 scripts/prepare-rdf.py"
    exit 1
fi

echo "Starting Apache Jena Fuseki..."
echo "  Data: $DATA_FILE"
echo "  Binding: localhost only (secure)"

cd "$FUSEKI_HOME"
java -Xmx4g -jar fuseki-server.jar \
    --localhost \
    --file="$DATA_FILE" \
    /gene > "$LOG_FILE" 2>&1 &

PID=$!
echo $PID > "$PID_FILE"

# Wait for startup
sleep 5

# Check if running
if kill -0 "$PID" 2>/dev/null; then
    TRIPLES=$(curl -s http://localhost:3030/gene/sparql \
        -H "Accept: application/json" \
        --data-urlencode 'query=SELECT (COUNT(*) AS ?c) WHERE { ?s ?p ?o }' \
        | grep -o '"value":"[0-9]*"' | head -1 | grep -o '[0-9]*')

    echo ""
    echo "Fuseki started successfully!"
    echo "  PID: $PID"
    echo "  Triples: ${TRIPLES:-unknown}"
    echo "  Endpoint: http://localhost:3030/gene/sparql"
    echo "  UI: http://localhost:3030 (localhost only)"
    echo ""
    echo "Query example:"
    echo "  ./scripts/sparql-query.sh 'SELECT * WHERE { ?s ?p ?o } LIMIT 5'"
    echo ""
    echo "Stop server:"
    echo "  ./scripts/stop-fuseki.sh"
else
    echo "Failed to start Fuseki. Check log: $LOG_FILE"
    cat "$LOG_FILE" | tail -20
    exit 1
fi
